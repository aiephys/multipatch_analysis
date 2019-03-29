"""Simple tool for displaying DB connections, used for debugging performance issues
"""
from __future__ import print_function, division
import os, sys, subprocess, re
from multipatch_analysis import config
from sqlalchemy import create_engine


engine = create_engine(config.synphys_db_host_rw + '/postgres?application_name=db_monitor')
conn = engine.connect()


def list_db_connections():
    tr = conn.begin()
    query = conn.execute("select client_addr, client_port, datname, query, state, usename, application_name, pid from pg_stat_activity;")
    result = query.fetchall()
    tr.rollback()
    return result


_known_hostnames = {}
def hostname(ip):
    global _known_hostnames
    if ip in _known_hostnames:
        return _known_hostnames[ip]
    try:
        host = subprocess.check_output(['host', ip]).partition('pointer ')[2].rstrip('.\n')
    except subprocess.CalledProcessError:
        host = "hostname not found"
    _known_hostnames[ip] = host
    return host


def check():
    print("====================  DB connections ======================")    
    connects = list_db_connections()
    connect_ips = [conn[0] for conn in connects]
    
    ips = list(set(connect_ips))
    counts = {ip:connect_ips.count(ip) for ip in ips}
    ips.sort(key=lambda ip: counts[ip], reverse=True)
    known_addrs = config.known_addrs
    
    for ip in ips:
        if ip is None:
            host = "[None]"
        else:
            host = hostname(ip)
        name = known_addrs.get(ip, known_addrs.get(host, '???'))
        count = counts[ip]
        print("{:10s}{:15s}".format(name, ip))
        
        for con in connects:
            if con[0] == ip:
                app = con.application_name
                for pkg in ['acq4', 'multipatch_analysis']:
                    a,b,c = app.partition(pkg)
                    if b == '':
                        continue
                    else:
                        app = c
                        break
                        
                state = '[%s]' % con.state
                query = con.query.replace('\n', ' ')[:120]
                
                print("          {:15s} {:15s} {:45s} {:6d} {:10s} {:s}   ".format(con.usename, con.datname, app[:45], con.pid, state, query))
    

class ScreenPrint(object):
    """Simple interface for printing full-screen text
    """
    def __init__(self):
        global print
        self._print = print
        print = self.print_line
        
        self.current_row = 0
        self.rows = 100
        self.columns = 100

    def reset(self):
        self.current_row = 0
        self._print("\033[0;0f")
        self.rows, self.columns = list(map(int, os.popen('stty size', 'r').read().split()))

    def print_line(self, msg):
        while len(msg) > self.columns:
            line, msg = msg[:self.columns], msg[self.columns:]
            self._print(line)
            self.current_row += 1
        self._print(("{:%ds}" % self.columns).format(msg))
        self.current_row += 1
        
    def clear_to_bottom(self):
        for i in range(self.current_row, self.rows-1):
            self.print_line("")


if __name__ == '__main__':
    import time
    
    sp = ScreenPrint()
    
    while True:
        sp.reset()
        check()
        sp.clear_to_bottom()        
        time.sleep(5)
