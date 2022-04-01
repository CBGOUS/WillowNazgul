# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:39:57 2021

@author: norab
"""

import sys
import yaml
from treeRun import treeRun
import timeit
import logging
import sys
import os
from datetime import datetime
import hashlib
import shortuuid


def initLogger(configFilepath):
    """
    initiate loggerfile
    :param configFilepath: path to configuration file
    :return:
    """

    # setup log folder based on configfile - folder + date + unique id
    logFolder = configFilepath.rsplit('/', 1)[0] + '/'

    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)

    now = datetime.now()
    dt_string = now.strftime("%d%m%Y_%H%M%S")
    id = shortuuid.ShortUUID().random(length=5)

    logfileName = os.path.join(logFolder, "logfile_" + dt_string + '_' + id + ".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)

    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)  # set the new handler
    log.addHandler(handler)
    logging.info("+" + "*" * 78 + "+")
    logging.info("project log file is <" + logfileName + ">")
    logging.info("+" + "*" * 78 + "+")
    logging.debug("debug mode is on")


def main(configFilepath=None):
    """
    Function: 
        - call function to run program
        - path to configuration file given as 1st argument
        
    """

    print("Job initiated")
    md5String = hashlib.md5(b"Nazgul").hexdigest()
    # Read configfile argument
    if configFilepath:
        configFilepath = configFilepath.strip()
    else:
        configFilepath = sys.argv[1]  # uncomment this to run from terminal

    with open(configFilepath, 'r') as c:
        configFile = yaml.safe_load(c)

    # initiate logger, logfile directed to output folder
    logFilepath = configFile.get('output_folder').strip()
    initLogger(logFilepath)
    # Run
    treeRun(configFile)


if __name__ == '__main__':
    pass
    #### Bash run
    # main()

    #### Run in IDE: 
    # Test job: 
    # configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/test_GDR_10genes/job_input/main_config.yml'

    # GDR calculation job:
    # configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/GDRcalculation/job_input/main_config.yml'

    # nullGDR calculation job:
    # configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/GDRcalculation/job_input/main_config.yml'

    # non-zero phydists job: 
    # configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/GDRcalculation/job_input/main_config.yml'

    # Run
    # starttime = timeit.default_timer()
    # main(configFilepath)
    # endtime = timeit.default_timer()

    # Testrun 25.03.22
    # configFilepath = 'C:/Users/norab/Master/WillowProject/testRuns/testRun25.03/main_config.yml'
    # main(configFilepath)

    # Testrun GDRnull 01.04.2022
