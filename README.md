Ns3-rmcat for GCC
=================

GCC(Google Congestion Control)
---------------------------------------------------------------
Algorithm is referenced in `draft-alvestrand-rmcat-congestion`


Implementation
--------------

To implement GCC in NS3 we refer to the WebRTC code in chromium. 


To clone WebRTC GCC repo:

    git clone https://chromium.googlesource.com/external/webrtc
     

Commit Hash ID of WebRTC's git which we refer is `3613fef7a2f0cae1b1e6352b690979d430626056`.

GCC module defines ns3 application (see `model/apps/gcc-node.cc <model/apps/gcc-node.cc>`_) running on network topologies. These ns3 application send fake video codec data (model/syncodecs, more on syncodecs), according to the GCC algorithm. Since `gcc-node` includes both receiver and sender functions, implementing half-duplex or full-duplex communication all is possible.

In half-duplex communication, sender node(gcc-node) sends fake video codec data in media packets to receiver node(gcc-node). Receiver node gets the sequence of packets and run delay-based congestion control algorithm(model/congestion-control/gcc-receiver-controller.cc). Receiver node sends a RTCP feedback packets for the loss information and a REMB feedback packets for the delay-based send bitrate. Sender node gets feedback packets from recevier, and run sender-based congestion control algorithm(model/congestion-control/gcc-sender-controller.cc) to get bandwidth estimation. Sender node then uses this bandwidth estimation to control the fake video encoder by adjusting its target video bitrate.

In full-duplex communication, each node(gcc-node) can act as both sender and receiver.


GCC Log Parser
--------------
To be announced.










Simulation
==========

Example 1 : rmcat-half-duplex-eval  
---------------------------------

### Example 1's Topology


```
   +---+                                                           +---+
   |S1 |====== \                 Forward -->              / =======|R1 |
   +---+       \\                                        //        +---+
                \\                                      //
   +---+       +-----+                               +-----+       +---+
   |S2 |=======|  A  |------------------------------>|  B  |=======|R2 |
   +---+       |     |<------------------------------|     |       +---+
               +-----+                               +-----+
   (...)         //                                     \\         (...)
                //          <-- Backward                 \\
   +---+       //                                         \\       +---+
   |Sn |====== /                                           \ ======|Rn |
   +---+                                                           +---+
```

### Default run
`./waf --run "rmcat-half-duplex-eval" > result.out 2> &1`
+ default
  + number of rmcat flows : 1
  + rmcat alogirhtm : CCFS
  + number of tcp : 0
  + number of udp : 0
  + throughput : 1000 Kbps
  + rtt : 300ms
  + qdelay : 300ms
  
  
### Related options

+ `--algo=nada/ccfs/gcc` : kind of rmcat flow (NADA, CCFS, GCC)
+ `--rmcat=N` : number of rmcat flows
+ `--tcp=N`: tcp cross traffic
+ `-udp=N`: udp cross traffic
+ `--numSession=N` : number of sessions per rmcat flow
+ `--kbps=xxxx` : config data rate option(Throuput)

### Run with options
`./waf --run "rmcat-half-duplex-eval --algo=ccfs --rmcat=1 --tcp=1"`

`./waf --run "rmcat-half-duplex-eval --algo=gcc --rmcat=3"`

### Result



Example 2 : rmcat-full-duplex-eval
----------------------------------
To be announced.



References
==========
1. draft-alvestrand-rmcat-congestion : (https://tools.ietf.org/html/draft-ietf-rmcat-gcc-02)
2. draft-sarker-rmcat-eval-test : (https://tools.ietf.org/html/draft-ietf-rmcat-eval-test-06)
3. 






