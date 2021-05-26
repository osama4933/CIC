# CIC-Code

Instructions to Use CIC demodulation on LoRa Packets

1- open param_configs file and set SF, BW, Receiver Sampling Frequency 
and num of data symbols in a LoRa pkt accroding to your LoRa pkts' and 
receiver configuration. set appropriate path to the binary file that 
contains I/Q samples of LoRa signal in the variable named 'path' & in the 
variable named 'fil_nm' set the appropriate binary file name  to be read 
by matlab.

2- Open main.m file
   -> main.m file has been divided into following sections
	- %%  Loading variables
	- %%  Loading File
	- %%  Active Sessions Detection using Dechirping Windows
	- %%  DC correlations to find LoRa pkts out of collision
	- %%  UC correlation to filter false positives and frequency 
		offset & Packets' SNR estimation
	- %%  Data Demodulation using CIC
	- %%  Calculating SER for the File

3- Run whole main file, or run the file section by section avoiding for loop iterating 
over 'm' (choose an appropriate active session yourself) (while 
running each section, comment 'continue' statements in 'if' conditions 
in sections '%%  DC correlations to find LoRa pkts out of collision' 
and '%%  UC correlation to filter false positives and frequency offset
 & Packets' SNR estimation' because 'continue' statement can only be called in a for loop)

4- In the current folder, 'sym.mat' is supposed to contain 'sym' variable 
that contains true symbols against which output of CIC is to be compared in 
order to compute 'Symbol Error Rate' or 'Symbol Throughput'
 ('symbol throughput calculation not implemnted')

5- The Symbol Error Rate is computed and printed on command window. At the end 
'demod_sym_stack' (2D matrix) variable can be found. The number of rows indicate 
the packets detected in the file. At each row, lies a demodulated packet whose 
column contain corresponding demodulated symbols.

6- the section named '%%  Active Sessions Detection ...' finds the sessions that 
have active LoRa tranmissions of a specific SF going on. If the Binary file has 
data of high network traffic, this may result in very long active sessions. Passing 
such sessions to DC_Correlation may take a lot of time. Therefore such long sessions 
should be chunked into smaller units with appropriate sample overlaps in between smaller 
chunks, if the length of any session exceeds some threshold length.


******************************DATA***************************************************
2 DataSets have been provided and can be accessed via following links
Low SNR Data: https://uwprod-my.sharepoint.com/:u:/g/personal/mshahid2_wisc_edu/EQlJgjRMH89Hp3kIEBkhwA8BsR_b2dqOQzNrjurTMyZI2g?e=A2MHJJ
High SNR Data: https://uwprod-my.sharepoint.com/:u:/g/personal/mshahid2_wisc_edu/EaSovqXfd8dEs3fy65DGqbMBEbpAyZz23Y-GrLGG73R_9Q?e=iVmbpa
Naming Convention: '# of LoRa Tx in network'_'individual node's pkt rate'lm_1
 e.g. 20tx_1lm_1 -> 20 = '# of LoRa Tx in network'
		    1lm = 'individual node's pkt rate'
Experimental Setup has been Described in the Experimental_setup.pdf.
