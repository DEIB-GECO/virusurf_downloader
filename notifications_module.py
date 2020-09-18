from loguru import logger
from notifiers.logging import NotificationHandler
import os

"""
Created by tomalf2 on Sep, 2020.

This module enables additional handlers for get someone notified in case of unexpected error conditions. 
This function makes use of third-parthy services requiring private authentication information. Authentication
parameters are stored in a separate file called "notifier_params.tsv" to protect the user of the system.

This module implements the configuration for using Telegram service only, but other sevices can be used too.
Refer to https://loguru.readthedocs.io/en/0.1.0/api/notifier.html for further details.   
"""


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


def setup_telegram_notifier():
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
        logger.add(notification_handler, level="ERROR")


def setup_any_additional_error_notifiers():
    try:
        setup_telegram_notifier()
    except:
        logger.warning("Telegram notifier not configured.")
