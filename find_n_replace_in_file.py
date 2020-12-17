"""
Created by tomalf2 on dic, 2020.
"""
import sys
import os

path_to_file = sys.argv[1]
word_to_find = sys.argv[2]
word_to_replace_with = sys.argv[3]

path_temp_file = path_to_file+".temp"

with open(path_to_file, mode="r") as original:
    with open(path_temp_file, mode="w") as modified:
        for line in original.readlines():
            new_line = line.replace(word_to_find, word_to_replace_with)
            modified.write(new_line)

os.remove(path_to_file)
os.rename(path_temp_file, path_to_file)