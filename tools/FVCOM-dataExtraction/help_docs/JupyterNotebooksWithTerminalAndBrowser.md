### Setting up the remote notebook
1.	In the terminal, SSH into a specific node on kestrel to host the session (e.g., kd3)
2.	Load the environment that you would like to use for the jupyter notebook (e.g., `module load conda && conda activate <env_name>`)
3.	Make sure you are in a location on kestrel that includes the notebook you wish to open.
4.	Launch the notebook capability using `jupyter notebook --no-browser --ip=*`
5.	The screen will print out a bunch of lines, showing at the bottom “Or copy and paste one of these URLs”. You will need to copy this after the following step

### Setting up the remote connection
6.	In a new terminal tab or window, open a connection to the same node with a dedicated port (e.g., `ssh -L 8888:kd3:8888 kestrel.hpc.nlr.gov`) 

### Opening the notebook
7.	Copy one of the URLs from step 5 and paste it into a browser. Sometimes only one of the two URLs works, so be sure to try the other if the first one does not open properly. 
8.	Navigate to and open your desired notebook on kestrel.
