## Project 0 Writeup
#### Bohao Tang

#### 1. How large is the system's root drive, mounted on "/"?
Answer: 8065444 kb
Command: `df /`

#### 2. How much memory does the system have?
Answer: 1014564 kb
Command: `free -k`

#### 3. How much memory is being used?
Answer: 57928 kb
Command: `free -k`

#### 4. What linux kernel version are you running?
Answer: 4.4.0-1047-aws
Command: `uname -r`

#### 5. What is the MAC address of the Ethernet card?
Answer: 12:d4:dc:af:cc:94
Command: `ifconfig -a`

#### 6. What is the last message in the log file where system messages/errors would be found? Where is this file located?
Answer: "Feb  4 18:17:01 ip-10-0-144-179 CRON[18263]: (root) CMD (   cd / && run-parts --report /etc/cron.hourly)" in file "syslog" in "/var/log"
Command: `tail /var/log/syslog`
