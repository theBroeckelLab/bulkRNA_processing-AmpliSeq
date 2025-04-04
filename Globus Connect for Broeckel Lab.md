Author: Ryan Gallagher (04-04-2025)

## Globus Connect for Broeckel Lab

This document will step the user through setting up Globus Connect. This software is used to transfer data from NWGS for our AmpliSeq results. This document will mainly follow the [RCC documentation](https://docs.rcc.mcw.edu/storage/globus/?h=globus#setup).

### 1. Setup
#### 1.1 Create a Key and Start Globus
Open your browser to being creating a [Collection](https://app.globus.org/collections/gcp?generate_key). On this page, create a collection named "NWGS Directories". Then, obtain your setup key. 

Now, login to the RCC and type the following in the command line:

`module load globusconnect`

`globusconnect -setup KEY_HERE`

Where `KEY_HERE` is your setup key. Upon success, enter the following into the RCC command line:

`globusconnect -start`

#### 1.2 Map Relevant Directories
Once we've configured our collection and started Globus, we will end the session and edit some paths.

Press `<Control+C>` or `<Ctrl + C>` to escape the Globus Session.

Next, in the RCC command line, navigate via:

`cd ~/.globusonline/lta/` 

And open the `config-paths` file:

`vim config-paths`

Within this file are the paths mapped to your collection. We will add the path:

```
~/,0,1
/scratch/g/ubroecke,0,1
/group/ubroecke,0,1
```

**Be sure you've spelled everything correctly!**



You can exit vim via `:wq`.

Go ahead and start Globus Connect again:

`globusconnect -start`

### 2. Connect to Globus Personal

While `globusconnect -start` is running, navigate to https://app.globus.org/file-manager or use the link provided by the lab from the NWGC files. 

