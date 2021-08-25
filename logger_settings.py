import logging
from typing import Optional

from loguru import logger
from notifiers.logging import NotificationHandler
import notifiers.exceptions
import os
from tqdm import tqdm
import warnings
import sys
import requests
import threading

"""
Created by tomalf2 on Sep, 2020.

This module enables additional handlers to get someone notified in case of unexpected error conditions. It can be used
also to explicitly send a notification whenever needed.
This function makes use of third-parthy services requiring private authentication information. Authentication
parameters are stored in a separate file called "notifier_params.tsv" to protect the user of the system.

This module implements the configuration for using Telegram service only, but other sevices can be used too.
Refer to https://loguru.readthedocs.io/en/0.1.0/api/notifier.html for further details.   
"""
messenger_parameters = dict()   # container for auth params of any messaging services


def _read_params():
    parameter_file_path = f".{os.path.sep}notifier_params.tsv"
    if not os.path.exists(parameter_file_path):
        with open(parameter_file_path, mode="w") as file:
            file.write(
                "# Lines starting with # are comments.\n"
                "# Example parameters:\n"
                "# telegram\t<API-TOKEN>\t<target_chat_id>\n"
                "# Every entry is separated by a TAB"
            )
            raise ValueError("You need to specify the parameters for receiving notifications in the file "
                             "./notifier_params.tsv in order to use this feature.")
    else:
        param_list = list()
        with open(parameter_file_path, mode="r") as file:
            for line in file:
                sline = line.rstrip().lstrip()
                if sline.startswith("#"):
                    continue
                param_list.append(sline.split())
        if not param_list:
            raise ValueError("You need to specify the parameters for receiving notifications in the file "
                             "./notifier_params.tsv in order to use this feature.")
        else:
            return param_list


def _setup_telegram_notifier():
    t_settings = [x for x in _read_params() if x[0].lower() == "telegram"]
    if not t_settings:
        raise ValueError("To use the telegram notifier you need to edit the file ./notifier_params.tsv and specify "
                         "'telegram' followed by your API_TOKEN and the target chat_id")
    for config in t_settings:
        notification_handler = NotificationHandler(
            provider="telegram", defaults={
                "token": config[1],
                "chat_id": config[2]
            }
        )
        notification_handler_wrapper = NotificationHandlerWrapper(notification_handler)
        logger.add(sink=notification_handler_wrapper,
                   level="ERROR",
                   format="<green>{time:YYYY-MM-DD HH:mm:ss Z}</green> | "
                          "<level>{level: <8}</level> | "
                          "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
                   backtrace=False, diagnose=False)  # disables sending exception stacktrace


def setup_any_additional_error_notifiers():
    try:
        _setup_telegram_notifier()
    except:
        logger.warning("Telegram notifier not configured.")


class NotificationHandlerWrapper(logging.Handler):
    """
    Wraps the original NotificationHandler from notifiers.logging and ensures that each message reaches at most the
    size of 4096 chars. Messages above such threshold are split into chunks and sent as separate messages
    """
    MAX_MSG_LENGTH = 4096
    MAX_MSG_LENGTH_WITH_TRAILING_ADDITION = 4096 -24

    def __init__(self, true_notification_handler):
        super().__init__()
        self.handler = true_notification_handler

    def emit(self, record: logging.LogRecord):
        # send up to 4096 chars per notification
        message: str = record.getMessage()
        try:
            while len(message) > NotificationHandlerWrapper.MAX_MSG_LENGTH:
                record.msg = f"{message[:NotificationHandlerWrapper.MAX_MSG_LENGTH_WITH_TRAILING_ADDITION]} -- CONTINUE IN NEXT MSG"
                self.handler.emit(record)
                message = message[NotificationHandlerWrapper.MAX_MSG_LENGTH_WITH_TRAILING_ADDITION:]
            record.msg = message
            self.handler.emit(record)
        except notifiers.exceptions.NotificationError:
            logger.info(f"Last attempt to send a notification message has been rejected because the sending rate was too high.")

    def handleError(self, record: logging.LogRecord):
        self.handler.handleError(record)

    def init_providers(self, provider, kwargs):
        self.handler.init_providers(provider, kwargs)


def setup_logger(log_file_prefix: str):
    logger.remove()  # removes default logger to stderr with level DEBUG
    # on console print from level INFO on
    logger.add(sink=lambda msg: tqdm.write(msg, end=''),
               level='TRACE',
               format="<green>{time:YYYY-MM-DD HH:mm:ss Z}</green> | "
                      "<level>{level: <8}</level> | "
                      "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
               colorize=True,
               backtrace=True,
               diagnose=True,
               enqueue=True)
    # log to file any message of any security level
    logger.add("./logs/" + log_file_prefix + "_{time}.log",
               level='TRACE',
               rotation='100 MB',
               format="<green>{time:YYYY-MM-DD HH:mm:ss Z}</green> | "
                      "<level>{level: <8}</level> | "
                      "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
               colorize=False,
               backtrace=True,
               diagnose=True,
               enqueue=True)
    setup_any_additional_error_notifiers()  # attach a sink to loguru.logger sinks
    _setup_messengers()                      # setup an external communication services to be invoked with send_message

    # redirect warnings
    def customwarn(message, category, filename, lineno, file=None, line=None):
        logger.warning(warnings.formatwarning(message, category, filename, lineno))

    warnings.showwarning = customwarn

    logger.info(f"EXECUTING main.py {' '.join(sys.argv[1:])} on PID {os.getpid()}")


def _setup_messengers():
    try:
        _setup_telegram_messenger()
    except:
        logger.warning("Telegram messenger not configured.")


def _setup_telegram_messenger():
    global messenger_parameters
    t_settings = [x for x in _read_params() if x[0].lower() == "telegram_messenger"]
    if not t_settings:
        raise ValueError("To use the telegram notifier you need to edit the file ./notifier_params.tsv and specify "
                         "'telegram' followed by your API_TOKEN and the target chat_id")
    else:
        telegram_messengers = messenger_parameters.get("telegram_messenger", [])
    for config in t_settings:
        param = {
            "token": config[1],
            "chat_id": config[2]
        }
        telegram_messengers.append(param)
    messenger_parameters["telegram_messenger"] = telegram_messengers


def _send_telegram_mesage(message: str):
    """
    Helper function to send a message using the available telegram_messenger configurations. Don't use this function
    directly; use send_message instead.
    """
    my_messengers = messenger_parameters.get("telegram_messenger", [])
    for a_messenger in my_messengers:
        try:
            requests.post(url=f"https://api.telegram.org/bot{a_messenger['token']}/sendMessage", timeout=5, data={
                "chat_id": a_messenger["chat_id"],
                "text": message
            })
        except requests.exceptions.Timeout:
            logger.error(f"Messenger of class telegram with token {a_messenger['token']} and chat_id "
                         f"{a_messenger['chat_id']} raised a timeout.")
        except:
            logger.exception(f"Error encounter when trying to use the telegram messenger with token {a_messenger['token']} "
                             f"and chat_id {a_messenger['chat_id']}.")


def send_message(message: str, block: bool = False) -> Optional[threading.Thread]:
    """
    Sends a message to the configured messengers (currently supports only telegram).
    If block is True, the network call will be done inside the thread calling this function, otherwise a new thread is
    spawn to take care of it and its object is returned (and can be joined if necessary). Nevertheless, the network call
    will either success or be cancelled after a timeout of 5 seconds.
    """
    if not block:
        thread = threading.Thread(target=_send_telegram_mesage, args=(message,))
        thread.start()
        return thread
    else:
        _send_telegram_mesage(message)
