#!/usr/bin/env python

import subprocess
import os
import re
import datetime
import csv
import multiprocessing
import time


def run(cmd):
    """Execute the specified command in a shell environment
    """
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         shell=True)

    # Wait for process completion
    (stdout, stderr) = p.communicate()
    retcode = p.returncode
    return stdout, retcode


def get_timestamp():
    """Return a universal timestamp for measurements."""
    timestamp = datetime.datetime.now()
    return timestamp


def take_snapshot():
    """Create dictionary of snapshots in the following format.
    FORMAT: {timestamp:9am, mic0_watts:300, mic1_watts:30, mic0_temp:89, 
    mic1_temp:98}
    """
    snapshot = {'timestamp': get_timestamp()}
    
    #Contact card for frequency and power data
    (power_out, retcode1) = run("micsmc -f")
    #Add the data for power and frequency to the dictionary (snapshot)
    if retcode1 == 0:
        snapshot.update(parse_micsmc("Total_Power", power_out, 
                                     "Total Power: ............. ",
                                     " Watts"))
        snapshot.update(parse_micsmc("Core_Freq", power_out, 
                                     "Core Frequency: .......... ", " GHz"))
    #Contact card for thermal data 
    (temp_out, retcode2) = run("micsmc -t")
    #Add the thermal data to the dictionary (snapshot)
    if retcode2 == 0:
        snapshot.update(parse_micsmc("CPU_Temp", temp_out, 
                                     "Cpu Temp: ................ ", " C"))
    
    return snapshot


def parse_micsmc(name, buff, before_key, after_key):
    """Parse the data returned from micsmc tool containing the thermal
    and power readings.
    """
    ret_dict = {}
    regex = "%s.*%s" % (re.escape(before_key), re.escape(after_key))
    mic = 0
    #add data to a dictionary to be returned
    for line in re.findall(regex, buff):
        re_obj = re.search("\d+\.\d+", line)
        ret_dict['%s_mic%d' % (name, mic)] = re_obj.group()
        mic += 1
    return ret_dict


def output_snapshots(filename, snapshots):
    """Write the given cummulative snapshot list to the filename given.
    snapshots - list of dictionaries [{}, {}, {}]
    """
    #Obtain dictionary keys and open csv file 
    fieldnames = snapshots[0].keys()
    with open(filename, "wb") as f:
        dw = csv.DictWriter(f, delimiter=",", fieldnames=fieldnames, 
                            quoting=csv.QUOTE_ALL)
        headers = {}

        #Add headers to csv file
        for n in dw.fieldnames:
            headers[n] = n
        dw.writerow(headers)

        #Add data to csv file
        for row in snapshots:
            dw.writerow(row)
    return


class Stats(multiprocessing.Process):
    """This class is used to control the process which is doing the
    logging of thermals and power of the mic cards.
    """
    def __init__(self, stat_file, snap_delay):
        multiprocessing.Process.__init__(self)
        self.exit = multiprocessing.Event()
        self.stat_file = stat_file
        self.snap_delay = snap_delay
        return

    def run(self):
        snap_list = []
        snap_time = 0
        while not self.exit.is_set():
            cur_time = time.time()
            if (cur_time - snap_time) > self.snap_delay:
                snap_time = cur_time
                snap_list.append(take_snapshot())
        output_snapshots(self.stat_file, snap_list)
        return

    def shutdown(self):
        self.exit.set()
        return


if __name__ == "__main__":
    f = "output.csv"
    stats = Stats(f, 0.5)

    stats.start()
    raw_input("Press enter to stop")

    stats.shutdown()
    stats.join()
