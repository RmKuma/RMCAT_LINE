/******************************************************************************
 * Copyright 2016-2017 cisco Systems, Inc.                                    *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License");            *
 * you may not use this file except in compliance with the License.           *
 * You may obtain a copy of the License at                                    *
 *                                                                            *
 *     http://www.apache.org/licenses/LICENSE-2.0                             *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 ******************************************************************************/

/**
 * @file
 * Simple example demonstrating the usage of the rmcat ns3 module, using:
 *  - NADA as controller for rmcat flows
 *  - Statistics-based traffic source as codec
 *  - [Optionally] TCP flows
 *  - [Optionally] UDP flows
 *
 * @version 0.1.1
 * @author Jiantao Fu
 * @author Sergio Mena
 * @author Xiaoqing Zhu
 */

#include "ns3/nada-controller.h"
#include "ns3/rmcat-sender.h"
#include "ns3/rmcat-receiver.h"
#include "ns3/rmcat-constants.h"
#include "ns3/point-to-point-helper.h"
#include "ns3/data-rate.h"
#include "ns3/bulk-send-helper.h"
#include "ns3/packet-sink-helper.h"
#include "ns3/udp-client-server-helper.h"
#include "ns3/internet-stack-helper.h"
#include "ns3/traffic-control-helper.h"
#include "ns3/ipv4-address-helper.h"
#include "ns3/core-module.h"
#include "ns3/ipv4-global-routing-helper.h"

#include "ns3/gcc-node.h"
#include <string>

// Maybe Ignore it 
const uint32_t GCC_DEFAULT_RMIN  =  150000;  // in bps: 150Kbps
const uint32_t GCC_DEFAULT_RMAX  = 1500000;  // in bps: 1.5Mbps

// TODO SHOULD MODIFY THIS BUT DON'T KNOW EXACT INITIAL VALUE.
const uint32_t GCC_DEFAULT_RINIT =  150000;  // in bps: 150Kbps (r_init)

const uint32_t TOPO_DEFAULT_BW     = 1000000;    // in bps: 1Mbps
const uint32_t TOPO_DEFAULT_PDELAY =      50;    // in ms:   50ms
const uint32_t TOPO_DEFAULT_QDELAY =     300;    // in ms:  300ms

using namespace ns3;

static void InstallTCP (Ptr<Node> sender,
                        Ptr<Node> receiver,
                        uint16_t port,
                        float startTime,
                        float stopTime)
{
    // configure TCP source/sender/client
    auto serverAddr = receiver->GetObject<Ipv4> ()->GetAddress (1,0).GetLocal ();
    BulkSendHelper source{"ns3::TcpSocketFactory",
                           InetSocketAddress{serverAddr, port}};
    // Set the amount of data to send in bytes. Zero is unlimited.
    source.SetAttribute ("MaxBytes", UintegerValue (0));
    source.SetAttribute ("SendSize", UintegerValue (DEFAULT_PACKET_SIZE));

    auto clientApps = source.Install (sender);
    clientApps.Start (Seconds (startTime));
    clientApps.Stop (Seconds (stopTime));

    // configure TCP sink/receiver/server
    PacketSinkHelper sink{"ns3::TcpSocketFactory",
                           InetSocketAddress{Ipv4Address::GetAny (), port}};

    auto serverApps = sink.Install (receiver);
    serverApps.Start (Seconds (startTime));
    serverApps.Stop (Seconds (stopTime));

}

static void InstallGccApps (Ptr<Node> node_1,
                         Ptr<Node> node_2,
                         uint16_t port_1,
                         uint16_t port_2,
                         float startTime,
                         float stopTime,
                         uint32_t numStream,
                         bool biDirec)
{
    Ptr<GccNode> app_1 = CreateObject<GccNode> ();
    Ptr<GccNode> app_2 = CreateObject<GccNode> ();
   
    node_1->AddApplication (app_1);
    node_2->AddApplication (app_2);
 
    app_1->SetUp (port_1,1000000); 
    
    if(biDirec)
      app_2->SetUp (port_2,1000000); 
    else
      app_2->SetUp (port_2,0); 
    
    Ptr<Ipv4> ipv4_1 = node_1->GetObject<Ipv4> ();
    Ptr<Ipv4> ipv4_2 = node_2->GetObject<Ipv4> ();

    Ipv4Address ipAdd_1 = ipv4_1->GetAddress (1, 0).GetLocal ();
    Ipv4Address ipAdd_2 = ipv4_2->GetAddress (1, 0).GetLocal ();

    app_1->SetDest(ipAdd_2, port_2);
    app_2->SetDest(ipAdd_1, port_1);

    uint64_t* streamSize = (uint64_t*)malloc(numStream*sizeof(uint64_t));

    for(uint32_t i=0;i<numStream;i++)
    {
      streamSize[i] = 1000000;
    }

    app_1->AddMulStream(numStream, streamSize);
    
    if(biDirec)
      app_2->AddMulStream(numStream, streamSize);
    else
      app_2->AddMulStream(0, NULL);
    
    const auto fps = 30.;		// Set Video Fps.
    auto innerCodec_1 = new syncodecs::StatisticsCodec{fps};
    auto innerCodec_2 = new syncodecs::StatisticsCodec{fps};
    auto codec_1 = new syncodecs::ShapedPacketizer{innerCodec_1, DEFAULT_PACKET_SIZE};
    auto codec_2 = new syncodecs::ShapedPacketizer{innerCodec_2, DEFAULT_PACKET_SIZE};
    app_1->SetCodec (std::shared_ptr<syncodecs::Codec>{codec_1});
    app_2->SetCodec (std::shared_ptr<syncodecs::Codec>{codec_2});

    app_1->SetStartTime (Seconds (startTime));
    app_1->SetStopTime (Seconds (stopTime));
    
    app_2->SetStartTime (Seconds (startTime));
    app_2->SetStopTime (Seconds (stopTime));
}

int main (int argc, char *argv[])
{
    // Number of Flows 
    int nRmcat = 1;
    std::string mode = "gcc";
    int nTcp = 0;
    int nUdp = 0;
    bool biDirec = false;
    float endTime = 1500;
    uint32_t numStream = 0;

    bool log = true;
    
    CommandLine cmd;
    cmd.AddValue ("rmcat", "Number of rmcat (GCC) flows", nRmcat);
    cmd.AddValue ("mode", "nada/gcc/vcc", mode);   // Default is declared in rmcat-sender.cc
    cmd.AddValue ("tcp", "Number of TCP flows", nTcp);
    cmd.AddValue ("udp",  "Number of UDP flows", nUdp);
    cmd.AddValue ("biDirec",  "Bi-Directional flow generation", biDirec);
    cmd.AddValue ("numStream",  "The numer of streams of Gcc nodes", numStream);
    cmd.AddValue ("log", "Turn on logs", log);
    cmd.Parse (argc, argv);

    if (log) {
        LogComponentEnable ("GccNode", LOG_LEVEL_ALL);
        LogComponentEnable ("GfpHeader", LOG_LEVEL_ALL);
        LogComponentEnable ("GccReceiverController", LOG_LEVEL_ALL);
        LogComponentEnable ("GccSenderController", LOG_LEVEL_ALL);
        LogComponentEnable ("TcpSocketBase", LOG_LEVEL_ALL);
    }

    // configure default TCP parameters
    Config::SetDefault ("ns3::TcpSocket::DelAckCount", UintegerValue (0));
    Config::SetDefault ("ns3::TcpL4Protocol::SocketType", StringValue ("ns3::TcpNewReno"));    // Tcp Type
    Config::SetDefault ("ns3::TcpSocket::SegmentSize", UintegerValue (1000));

    const uint64_t linkBw   = TOPO_DEFAULT_BW;
    const uint32_t msDelay  = TOPO_DEFAULT_PDELAY;
    const uint32_t msQdelay = TOPO_DEFAULT_QDELAY;

    NodeContainer nodes;
    nodes.Create(2);

    PointToPointHelper p2p;
    p2p.SetDeviceAttribute ("DataRate", DataRateValue  (DataRate (linkBw)));
    p2p.SetChannelAttribute ("Delay", TimeValue (MilliSeconds (msDelay)));
    auto bufSize = std::max<uint32_t> (DEFAULT_PACKET_SIZE, linkBw * msQdelay / 8000);
    p2p.SetQueue ("ns3::DropTailQueue",
                           "Mode", StringValue ("QUEUE_MODE_BYTES"),
                           "MaxBytes", UintegerValue (bufSize));
  
    NetDeviceContainer dev0 = p2p.Install (nodes);
    
    InternetStackHelper internet;
    internet.Install(nodes);
    
    Ipv4AddressHelper ipv4;
    ipv4.SetBase ("10.1.0.0", "255.255.255.0");
    ipv4.Assign (dev0);
    
    // disable tc for now, some bug in ns3 causes extra delay
    TrafficControlHelper tch;
    tch.Uninstall (dev0);
    
    int port = 8000;

    for (int i = 0; i < nRmcat; i++) {
        auto start = 100.*i;
        auto end = std::max (start + 1., endTime - start);

        uint16_t port_1 = port++;
        uint16_t port_2 = port++;

        NodeContainer endNodes;
        endNodes.Create(2);

        NodeContainer srcToRouter;
        srcToRouter.Add(endNodes.Get(0));
        srcToRouter.Add(nodes.Get(0));

        NodeContainer routerToRecv;
        routerToRecv.Add(nodes.Get(1));
        routerToRecv.Add(endNodes.Get(1));

        NetDeviceContainer srcDev = p2p.Install(srcToRouter);
        NetDeviceContainer recvDev = p2p.Install(routerToRecv);
        
        internet.Install(endNodes);

        std::string srcIpBase;
        srcIpBase = "10.1."+std::to_string(2*i+1)+".0";
        ipv4.SetBase(srcIpBase.c_str(), "255.255.255.0"); 
        ipv4.Assign(srcDev);

        std::string recvIpBase;
        recvIpBase = "10.1."+std::to_string(2*i+2)+".0";
        ipv4.SetBase(recvIpBase.c_str(), "255.255.255.0"); 
        ipv4.Assign(recvDev);
        
        tch.Uninstall(srcDev);
        tch.Uninstall(recvDev);
        
        if(mode == "gcc")
          InstallGccApps (endNodes.Get (0), endNodes.Get (1), port_1, port_2, start, end, numStream, biDirec);
    }

    for (int i = 0; i < nTcp; i++) {
        auto start = 101+ 10.*i;
        auto end = std::max (start + 1., endTime - start);
        uint16_t port_1 = port++;

        NodeContainer endNodes;
        endNodes.Create(2);

        NodeContainer srcToRouter;
        srcToRouter.Add(endNodes.Get(0));
        srcToRouter.Add(nodes.Get(0));

        NodeContainer routerToRecv;
        routerToRecv.Add(nodes.Get(1));
        routerToRecv.Add(endNodes.Get(1));

        NetDeviceContainer srcDev = p2p.Install(srcToRouter);
        NetDeviceContainer recvDev = p2p.Install(routerToRecv);
        
        internet.Install(endNodes);
        
        std::string srcIpBase;
        srcIpBase = "10.1."+std::to_string(2*(nRmcat+i)+1)+".0";
        ipv4.SetBase(srcIpBase.c_str(), "255.255.255.0"); 
        ipv4.Assign(srcDev);

        std::string recvIpBase;
        recvIpBase = "10.1."+std::to_string(2*(nRmcat+i)+2)+".0";
        ipv4.SetBase(recvIpBase.c_str(), "255.255.255.0"); 
        ipv4.Assign(recvDev);
        
        tch.Uninstall(srcDev);
        tch.Uninstall(recvDev);
        
        InstallTCP (endNodes.Get (0), endNodes.Get (1), port_1, start, end);

        if(biDirec)
        {
          uint16_t port_2 = port++;
          InstallTCP (endNodes.Get (1), endNodes.Get (0), port_2, start, end);
        }
    }
    
    Ipv4GlobalRoutingHelper::PopulateRoutingTables ();

    std::cout << "Running Simulation..." << std::endl;
    Simulator::Stop (Seconds (endTime));
    Simulator::Run ();
    Simulator::Destroy ();
    std::cout << "Done" << std::endl;

    return 0;
}
