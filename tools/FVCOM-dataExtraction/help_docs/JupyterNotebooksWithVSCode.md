### Setting up the host (Needed for the first time)
1.	Open VS Code
2.	If not already installed, install the extension “Remote – SSH” from Microsoft
3.	Go to the top (navigation bar) to type a command. Either type “>” or select “Show and Run Commands”
4.	Type (or click) the command “Remote-SSH: Connect to Host …”
5.	Select “+ Add New SSH Host”
6.	To add a new host, type the ssh command needed. For example, I want to run jupyter notebooks on kd3, so I will type `ssh <username>@kd3.hpc.nlr.gov`.
7.	After entering that, select the ssh configuration file you want to update. Typically you just want to update the one at `~/.ssh/config` (which would come up as `/Users/<username>/.ssh/config`).
8.	If you want, you can open the config file after saving it to give it a shorter label, like changing `Host kd3.hpc.nlr.gov` to `Host kd3`.

### Connecting to the host, opening and running jupyter notebook
1.	Open VS Code (if not already open)
2.	Go to the top (navigation bar) to type a command. Either type “>” or select “Show and Run Commands”
3.	Type (or click) the command “Remote-SSH: Connect to Host …”
4.	Select your desired host profile to connect (e.g., kd3). The window will connect, which might take a second. You must have an SSH key set up with the remote machine. If this is the first time you have connected to this remote machine from your local machine, you will need to affirm that you want to continue connecting.
5.	Now, in the new connected window that VS Code has opened, you will see SSH: at the bottom lest along with the host you are connected to. The next step is to open a folder on the remote machine. Click “Open Folder” in the left panel. Navigate to a folder containing the jupyter notebook you would like to run.
6.	Open the terminal (Cmd + \`), which will be on the remote machine. Do the steps needed to access the modules that the jupyter notebook needs (e.g., `module load conda && conda activate <env_name>`).
7.	Open the jupyter notebook from the terminal using `code /path/to/notebook.ipynb`.
8.	Go to the top right to “Select Kernel”. Pick “Jupyter Kernel”, and you should be able to find the python kernel associated with the active environment.
9.	Now the notebook is ready to run!
