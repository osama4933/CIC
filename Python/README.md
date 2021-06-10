# CIC Python
Instructions to use CIC demodulation on LoRa packets (Python version):

1. Run `pip install -r requirements.txt`.
1. Open `param_configs.py` and set SF, BW, Receiver Sampling Frequency and num of data symbols in a LoRa pkt accroding to your LoRa pkts' and receiver configuration. 
1. Open `main.py` and set appropriate path to the binary file that contains I/Q samples of LoRa signal in the variable named `in_file_path`.
1. (Optional) In `main.py`, set `symbols_ground_truth_path` to the path to the text file containing the expected symbols in the input file. The ground truth file should contain rows of space-separated symbols, where each row corresponds to a packet. See `symbols_ground.txt` for an example of the format.
1. Run `python main.py`. Demodulated symbols will be written to the file located at `symbols_out_path`, and peaks will be written to the file `peaks_out_path`.