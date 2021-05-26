# CIC-Matlab

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

5- The Symbol Error Rate/Symbol Throughput is computed and printed on command window. At the end, 
'demod_sym_stack' (2D matrix) variable can be found. The number of rows indicate 
the packets detected in the file. At each row, lies a demodulated packet whose 
column contain corresponding demodulated symbols.



******************************DATA***************************************************
3 DataSets have been provided and can be accessed via following link
https://uwprod-my.sharepoint.com/:f:/g/personal/mshahid2_wisc_edu/Es8DfWbahW1CpBdjv2uSdA8BdVTdkf0BblOOB472O8XtOQ

Naming Convention: '# of LoRa Tx in network'_'individual node's pkt rate'lm_1
 e.g. 20tx_1lm_1 -> 20 = '# of LoRa Tx in network'
		    1lm = 'individual node's pkt rate'
Experimental Setup has been Described in the Experimental_setup.pdf.







# CIC-Python

All the above instructions apply for Python version of CIC as well.



# Decoder

-> For Pkt throughput, we have used rpp0/gr-lora that can be installed along with dependencies as mentioned 
in the following 'rpp0' repository
https://github.com/rpp0/gr-lora

-> gr-lora contains two modules. 1- Standard LoRa Demodulator and 2- Decoder
LoRa Demodulator takes in I/Q samples and demodulates these samples based on default LoRa to give symbols
at the output, these symbols are then passed to Decoder that gives the decoded bits at the output. Since, We have 
our own CIC-based LoRa Demodulator, we run CIC to get the symbols and then pass our symbols output 
to gr-lora's decoder to get bits at the output (thus, bypassing gr-lora's demodulator).
The Code to bypass the gr-lora's demodulator can be found in the 'Decoder' folder.
Look for 2 files'decoder_impl_mod.cc' & 'decoder_impl_mod.h'.

-> After setting up the gr-lora move to folder /gr-lora/lib/. This contains 2 files 'decoder_impl.cc' and 'decoder_impl.h'. 
Copy the contents of 2 files in the Decoder folder of CIC and then paste it to the above mentioned files (copy contents 
of .cc to .cc and .h to .h files).

-> This modified version reads a text file as an input that contains demodulated symbols from CIC (file format: each row 
of txt file contains a single symbol and packets are separated by -1 identifier) and outputs decoded bits 
corresponding to each packet on a row in a text file.

-> Add correct paths to the input and output text files in the LoRa-Class constructor.

-> Get the Decoded bits in the relevant bits output file.