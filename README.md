# stealth-susy-dilepton
For development of dilepton channnel of Stealth/RPV SUSY search

To check this out, follow the basic git instructions.  For example, in a terminal window, you try typing this from the command line
```
git clone git@github.com:owen234/stealth-susy-dilepton.git
```
Once you have it checked out and [Root](https://root.cern) installed and working from the terminal, you can try using the code.  For example, start a root session from the command line and enter these lines:
```
root -l
root [0] .L dump1.c
root [1] dump1 d("root-files/signalR2-rpv_stop_750_t3j_uds.root")
root [2] d.Loop()
```
The "root \[\*\]" part in the lines above is the root prompt, not part of what you type.   It assumes that you put the root files in a subdirectory of the current directory named **root-files**.
That should open up three canvas windows on your screen and also print out a lot of information about the first event.

You can also try the other code like this
```
root -l
root [0] .L rpv_analysis2.c
root [1] rpv_analysis2 r("root-files/signalR2-rpv_stop_750_t3j_uds.root")
root [2] r.Loop()
root [3] h_rec_ht->Draw()
```

When Loop finishes, it prints the long list of histograms that were created.  The last command above draws one of them in a canvas.
