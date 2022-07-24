import os
from shutil import copy

# An XP folder contains always a log folder with all the logs
# A file named
# log-toto.txt
# in $(xp.path)/logs/
# will be copied in xp.dest folder
# and named :
# $(xp.name)-toto.txt


class XP:
    def __init__(self, path, name, destination):
        self.path = path  # Where the log folder in contained
        self.name = name  # Name to give as a prefix to the logs
        self.dest = destination  # Where to copy the files


def final_log_name(c_logname, prefix):
    return prefix + "-" + c_logname.split("-")[1]


# Where to find all the experimentations folders
global_path = "../Xp/"

# Where to copy all the log files
dest = "../XP/dest/"

# All the experimentations
XPs = [
    XP(global_path + "test/", "test", dest),
    XP(global_path + "prout/", "prout", dest)
]

for xp in XPs:
    for log in os.listdir(xp.path + "logs/"):
        copy(xp.path + "logs/" + log, xp.dest + final_log_name(log, xp.name))
