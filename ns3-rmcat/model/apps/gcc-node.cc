/******************************************************************************
 * Copyright 2016-2017 Cisco Systems, Inc.                                    *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License");            *
 * you may not use this file except in compliance with the License.           *
 *                                                                            *
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
 * Sender application implementation for rmcat ns3 module.
 *
 * @version 0.1.1
 * @author Jiantao Fu
 * @author Sergio Mena
 * @author Xiaoqing Zhu
 */

#include "gcc-node.h"
#include "ns3/udp-socket-factory.h"
#include "ns3/packet.h"
#include "ns3/simulator.h"
#include "ns3/uinteger.h"
#include "ns3/log.h"

#include <sys/stat.h>

NS_LOG_COMPONENT_DEFINE ("GccNode");
#define INITRATE 500000 //1Mbps
#define VIDEOINTERVAL 1000 //1s == 1000ms
#define AUDIOINTERVAL 5000 //5s == 5000ms
#define REMBINTERVAL 1000 //1s == 1000ms

#define THROCHECKINTERVAL 1000 //1s == 1000ms
#define DELAYCHECKINTERVAL 100 //100ms

namespace ns3 {

GccNode::GccNode ()
: m_destIp{}
, m_destPort{0}
, m_localPort{0}
, m_socket{NULL}
, m_localSsrc{rand()}
, m_numSrcSsrc{1}
, m_enqueueEvent{}
, m_sendEvent{}
, m_sendOversleepEvent{}
, m_rtcpEvent{}
, m_rembEvent{}
, m_rSend{INITRATE}
, m_totalRateShapingBuf{0}
, m_nextEnqSsrcIndex{0}
, m_nextSendSsrcIndex{0}
, m_sending{false}
, m_remoteSsrc{0}
, m_receiving{false}
, m_pDelay{0}
, m_recvDataBytes{0}
, m_recvAllBytes{0}
{}

GccNode::~GccNode () {}

void GccNode::SetCodec (std::shared_ptr<syncodecs::Codec> codec)
{
    m_codec = codec;
}

// TODO (deferred): allow flexible input of video traffic trace path via config file, etc.
void GccNode::SetCodecType (SyncodecType codecType)
{
    syncodecs::Codec* codec = NULL;
    switch (codecType) {
        case SYNCODEC_TYPE_PERFECT:
        {
            codec = new syncodecs::PerfectCodec{DEFAULT_PACKET_SIZE};
            break;
        }
        case SYNCODEC_TYPE_FIXFPS:
        {
            const auto fps = SYNCODEC_DEFAULT_FPS;
            auto innerCodec = new syncodecs::SimpleFpsBasedCodec{fps};
            codec = new syncodecs::ShapedPacketizer{innerCodec, DEFAULT_PACKET_SIZE};
            break;
        }
        case SYNCODEC_TYPE_STATS:
        {
            const auto fps = SYNCODEC_DEFAULT_FPS;
            auto innerStCodec = new syncodecs::StatisticsCodec{fps};
            codec = new syncodecs::ShapedPacketizer{innerStCodec, DEFAULT_PACKET_SIZE};
            break;
        }
        case SYNCODEC_TYPE_TRACE:
        case SYNCODEC_TYPE_HYBRID:
        {
            const std::vector<std::string> candidatePaths = {
                ".",      // If run from top directory (e.g., with gdb), from ns-3.26/
                "../",    // If run from with test_new.py with designated directory, from ns-3.26/2017-xyz/
                "../..",  // If run with test.py, from ns-3.26/testpy-output/201...
            };

            const std::string traceSubDir{"src/ns3-rmcat/model/syncodecs/video_traces/chat_firefox_h264"};
            std::string traceDir{};

            for (auto c : candidatePaths) {
                std::ostringstream currPathOss;
                currPathOss << c << "/" << traceSubDir;
                struct stat buffer;
                if (::stat (currPathOss.str ().c_str (), &buffer) == 0) {
                    //filename exists
                    traceDir = currPathOss.str ();
                    break;
                }
            }

            NS_ASSERT_MSG (!traceDir.empty (), "Traces file not found in candidate paths");

            auto filePrefix = "chat";
            auto innerCodec = (codecType == SYNCODEC_TYPE_TRACE) ?
                                 new syncodecs::TraceBasedCodecWithScaling{
                                    traceDir,        // path to traces directory
                                    filePrefix,      // video filename
                                    SYNCODEC_DEFAULT_FPS,             // Default FPS: 30fps
                                    true} :          // fixed mode: image resolution doesn't change
                                 new syncodecs::HybridCodec{
                                    traceDir,        // path to traces directory
                                    filePrefix,      // video filename
                                    SYNCODEC_DEFAULT_FPS,             // Default FPS: 30fps
                                    true};           // fixed mode: image resolution doesn't change

            codec = new syncodecs::ShapedPacketizer{innerCodec, DEFAULT_PACKET_SIZE};
            break;
        }
        case SYNCODEC_TYPE_SHARING:
        {
            auto innerShCodec = new syncodecs::SimpleContentSharingCodec{};
            codec = new syncodecs::ShapedPacketizer{innerShCodec, DEFAULT_PACKET_SIZE};
            break;
        }
        default:  // defaults to perfect codec
            codec = new syncodecs::PerfectCodec{DEFAULT_PACKET_SIZE};
    }

    // update member variable
    m_codec = std::shared_ptr<syncodecs::Codec>{codec};
}

void GccNode::SetDest(Ipv4Address dest_ip, uint16_t dest_port)
{
  m_destIp = dest_ip;
  m_destPort = dest_port;
}

void GccNode::SetUp (uint16_t local_port, uint64_t stream_size)
{
    if (!m_codec) {
        m_codec = std::make_shared<syncodecs::PerfectCodec> (DEFAULT_PACKET_SIZE);
    }

    m_recvController = std::make_shared<rmcat::GccRecvController> ();
    m_senderController = std::make_shared<rmcat::GccSenderController> ();

    m_localPort = local_port;
    m_socket = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());

    auto local = InetSocketAddress{Ipv4Address::GetAny (), local_port};
    auto ret = m_socket->Bind (local);
    NS_ASSERT (ret == 0);

    m_maxSize[m_localSsrc] = stream_size;

    m_socket->SetRecvCallback (MakeCallback (&GccNode::RecvPacket, this));
}

bool GccNode::AddMulStream(uint32_t num, uint64_t* stream_size)
{
  if(num > RembHeader::MAX_SSRCS_NUM)
    return false;

  for(uint32_t i=0;i<num;i++)
  {
    //uint32_t new_ssrc = rand(); //optional (not affect GCC performance)
    uint32_t new_ssrc = m_localSsrc + i;

    m_srcSsrcSet.insert(new_ssrc);
    m_maxSize[new_ssrc] = stream_size[i];
  }

  m_numSrcSsrc += m_srcSsrcSet.size();

  return true;
}

void GccNode::StartApplication ()
{
    // RTP initial values for sequence number and timestamp SHOULD be random (RFC 3550)
    m_sequence[m_localSsrc] = rand ();
    m_rtpTsOffset[m_localSsrc] = rand ();
    m_rateShapingBytes[m_localSsrc] = 0;
    m_enqBytes[m_localSsrc] = 0;
    m_complete = false;
    if(m_maxSize[m_localSsrc] == 0) //only receive mode 
      m_complete = true;


    for(const auto& ssrc : m_srcSsrcSet)
    {
      m_sequence[ssrc] = rand();
      m_rtpTsOffset[ssrc] = rand();
      m_rateShapingBytes[ssrc] = 0;
      m_enqBytes[ssrc] = 0;
    }

    m_rSend = INITRATE;

    m_nextEnqSsrcIndex = 0;
    m_nextSendSsrcIndex = 0;

    m_enqueueEvent = Simulator::Schedule (Seconds (0.0), &GccNode::EnqueuePacket, this);
    m_rtcpEvent = Simulator::Schedule(ns3::MilliSeconds(GetNextRtcpTime()), &GccNode::CreateRtcp, this);
    m_rembEvent = Simulator::Schedule (ns3::MilliSeconds(REMBINTERVAL), &GccNode::SendRemb, this); //REMB is sent per 1s (webRTC)

    uint64_t nowMs = Simulator::Now().GetMilliSeconds();
    m_lastThroCheckTime = ns3::MilliSeconds(uint64_t(nowMs/THROCHECKINTERVAL) * THROCHECKINTERVAL);
    m_lastDelayCheckTime = ns3::MilliSeconds(uint64_t(nowMs/DELAYCHECKINTERVAL) * DELAYCHECKINTERVAL);
    
}

void GccNode::StopApplication ()
{
  
    Simulator::Cancel (m_enqueueEvent);
    Simulator::Cancel (m_sendEvent);
    Simulator::Cancel (m_sendOversleepEvent);
    Simulator::Cancel (m_rtcpEvent);
    Simulator::Cancel (m_rembEvent);

    m_rateShapingBuf[m_localSsrc].clear ();
    m_rateShapingBytes[m_localSsrc] = 0;

    for(const auto& ssrc : m_srcSsrcSet)
    {
      m_rateShapingBuf[ssrc].clear();
      m_rateShapingBytes[ssrc] = 0;
    }
}

/*Send Methods*/
void GccNode::EnqueuePacket ()
{
    syncodecs::Codec& codec = *m_codec;
    codec.setTargetRate (m_rSend);
    ++codec; // Advance codec/packetizer to next frame/packet
    const auto bytesToSend = codec->first.size ();
    NS_ASSERT (bytesToSend > 0);
    NS_ASSERT (bytesToSend <= DEFAULT_PACKET_SIZE);

    uint32_t enq_ssrc = 0;
    bool canEnq = false;

    for(uint32_t i = 0;i<m_numSrcSsrc;i++)
    {
      if(m_nextEnqSsrcIndex == 0)
        enq_ssrc = m_localSsrc;
      else
      {
        std::set<uint32_t>::iterator it = m_srcSsrcSet.begin();
        for(uint32_t i=0;i<m_nextEnqSsrcIndex-1;i++)
          it++;
      
        enq_ssrc = *it;
      }

      if(m_enqBytes[enq_ssrc] >= m_maxSize[enq_ssrc])
      {
        m_nextEnqSsrcIndex++;
        if(m_nextEnqSsrcIndex == m_srcSsrcSet.size()+1)
          m_nextEnqSsrcIndex = 0;

        continue;
      }
      else
      {
        canEnq = true;
        break;
      }
          
    }

    if(!canEnq)
    {
      //do nothing...
      return;
    }

    m_rateShapingBuf[enq_ssrc].push_back (bytesToSend);
    m_rateShapingBytes[enq_ssrc] += bytesToSend;
    m_enqBytes[enq_ssrc] += bytesToSend;
    m_totalRateShapingBuf++;

    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::EnqueuePacket, packet enqueued, packet length: " << bytesToSend
                 << ", buffer size: " << m_rateShapingBuf[enq_ssrc].size ()
                 << ", buffer bytes: " << m_rateShapingBytes[enq_ssrc] << ", total enqueue bytes : "<<m_enqBytes[enq_ssrc]);

    double secsToNextEnqPacket = codec->second;
    Time tNext{Seconds (secsToNextEnqPacket)};
    m_enqueueEvent = Simulator::Schedule (tNext, &GccNode::EnqueuePacket, this);

    if (!m_sendEvent.IsRunning()) //send Start
    {
      m_sendEvent = Simulator::Schedule (tNext, &GccNode::SendPacket, this, 0);;
    }
}

void GccNode::SendPacket (uint64_t msSlept)
{
    uint32_t send_ssrc = 0;
    bool can_send = false;
    for(uint32_t i=0;i<m_numSrcSsrc;i++)
    {
      if(m_nextSendSsrcIndex == 0)
        send_ssrc = m_localSsrc;
      else
      {
        std::set<uint32_t>::iterator it = m_srcSsrcSet.begin();
        for(uint32_t i=0;i<m_nextEnqSsrcIndex-1;i++)
          it++;
        send_ssrc = *it;
      }

      if(m_rateShapingBuf[send_ssrc].size() == 0)
      {
        if(m_enqBytes[send_ssrc] >= m_maxSize[send_ssrc] && m_byeSet.find(send_ssrc) == m_byeSet.end())
          m_byeSet.insert(send_ssrc);
        
        m_nextSendSsrcIndex++;
        if(m_nextSendSsrcIndex == m_srcSsrcSet.size()+1)
          m_nextSendSsrcIndex = 0;
        continue;
      }
      else
      {
        can_send=true;
        break;
      }
    }

    if(!can_send)
    {
      NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"All stream buffer is empty.");
      return;
    }
   
    NS_ASSERT (m_rateShapingBuf[send_ssrc].size () > 0);
    NS_ASSERT (m_rateShapingBytes[send_ssrc] < MAX_QUEUE_SIZE_SANITY);

    const auto bytesToSend = m_rateShapingBuf[send_ssrc].front ();
    NS_ASSERT (bytesToSend > 0);
    NS_ASSERT (bytesToSend <= DEFAULT_PACKET_SIZE);
    m_rateShapingBuf[send_ssrc].pop_front ();
    m_totalRateShapingBuf--;
    NS_ASSERT (m_rateShapingBytes[send_ssrc] >= bytesToSend);
    m_rateShapingBytes[send_ssrc] -= bytesToSend;

    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::SendPacket, packet dequeued, packet length: " << bytesToSend
                 << ", buffer size: " << m_rateShapingBuf[send_ssrc].size ()
                 << ", buffer bytes: " << m_rateShapingBytes[send_ssrc]);

    // Synthetic oversleep: random uniform [0% .. 1%]
    uint64_t oversleepMs = msSlept * (rand () % 100) / 10000;
    Time tOver{MilliSeconds (oversleepMs)};
    m_sendOversleepEvent = Simulator::Schedule (tOver, &GccNode::SendOverSleep,
                                                this, bytesToSend,send_ssrc);

    // schedule next sendData
    const double msToNextSentPacketD = double (bytesToSend) * 8. * 1000. / m_rSend;
    const uint64_t msToNextSentPacket = uint64_t (msToNextSentPacketD);

    Time tNext{MilliSeconds (msToNextSentPacket)};
    m_sendEvent = Simulator::Schedule (tNext, &GccNode::SendPacket, this, msToNextSentPacket);
}

void GccNode::SendOverSleep (uint32_t bytesToSend, uint32_t send_ssrc) 
{
    if(!m_sending)
      m_sending = true;
    
    const auto nowMs = Simulator::Now ().GetMilliSeconds ();
    //m_sendController->... required

    ns3::RtpHeader header{96}; // 96: dynamic payload type, according to RFC 3551
    header.SetSequence (m_sequence[send_ssrc]++);
    NS_ASSERT (nowMs >= 0);
    // Most video payload types in RFC 3551, Table 5, use a 90 KHz clock
    // Therefore, assuming 90 KHz clock for RTP timestamps
    // header.SetTimestamp (m_rtpTsOffset[send_ssrc] + uint32_t (nowMs * 90 / 1000));;
    header.SetTimestamp (uint32_t(nowMs)); //we use simulator time as timestamp for easy analysis
    header.SetSsrc (send_ssrc);

    auto packet = Create<Packet> (bytesToSend);
    packet->AddHeader (header);

    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::SendOverSleep, " << packet->ToString ());
    m_socket->SendTo (packet, 0, InetSocketAddress{m_destIp, m_destPort});
}


void GccNode::CreateRtcp()
{
    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode:CreateRtcp");
    GccRtcpHeader* header;
    if(m_sending)
    {
      header = new GccRtcpHeader(GccRtcpHeader::RTCP_SR);
      header->SetSendSsrc(m_localSsrc);
      
      NS_ASSERT(header->AddSRFeedback(0,0,0,0,0) == GccRtcpHeader::RTCP_NONE); //Data in SR does not affect Gcc Performance
    }
    else
      header = new GccRtcpHeader(GccRtcpHeader::RTCP_RR);

    uint32_t rrSize = m_recvSsrcSet.size();

    for(const auto& ssrc : m_recvSsrcSet)
    {
      double fractionRatio = m_lost[ssrc]/(double)m_recvPackets[ssrc];
      uint8_t frac = uint8_t(255.0 * fractionRatio);
      uint32_t ret = 0;
      uint32_t nowMs = Simulator::Now().GetMilliSeconds();
      NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::Add RrBlock, Size "<<rrSize);

      if(m_lastSrRecvTime.find(ssrc) == m_lastSrRecvTime.end())
        m_lastSrRecvTime[ssrc] = 0;
      ret = header->AddRRFeedback(ssrc, frac, m_cumLost[ssrc], m_recvSeq[ssrc], 0, m_lastSrRecvTime[ssrc], nowMs-m_lastSrRecvTime[ssrc]);

      if(ret == GccRtcpHeader::RTCP_TOO_LONG)
      {
        header->SetTypeOrCount(31); //maximum number of RR blocks by RFC & WebRTC
        rrSize -= 31;
        SendRtcp(*header, false);
        
        delete header;
        header = new GccRtcpHeader(GccRtcpHeader::RTCP_RR);
        ret = header->AddRRFeedback(ssrc, frac, m_cumLost[ssrc], m_recvSeq[ssrc], 0, m_lastSrRecvTime[ssrc], nowMs-m_lastSrRecvTime[ssrc]);
      }
      NS_ASSERT(ret == GccRtcpHeader::RTCP_NONE);
      m_lost[ssrc] = 0;
      m_recvPackets[ssrc] = 0;
    }

    header->SetTypeOrCount(rrSize);

    if(!m_byeSet.empty())
    { 
      NS_ASSERT(header->AddByeFeedback(m_byeSet) == GccRtcpHeader::RTCP_NONE);
      NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"Add Bye Set, size = "<<m_byeSet.size());
      
      for(auto const& b_ssrc : m_byeSet)
      {
        if(b_ssrc != m_localSsrc)
        {
          m_srcSsrcSet.erase(b_ssrc);
          m_sequence.erase(b_ssrc);
          m_rtpTsOffset.erase(b_ssrc);
          m_rateShapingBuf.erase(b_ssrc);
          m_rateShapingBytes.erase(b_ssrc);
          m_enqBytes.erase(b_ssrc);
          m_maxSize.erase(b_ssrc);
          m_numSrcSsrc--;
        }
        else
        {
          m_sequence[b_ssrc] = 0;
          m_rtpTsOffset[b_ssrc] = 0;
          m_rateShapingBuf[b_ssrc].clear();
          m_rateShapingBytes[b_ssrc] = 0;
          m_enqBytes[b_ssrc] = 0;
          m_maxSize[b_ssrc] = 0;
          m_complete = true;
        }
      }
      m_byeSet.clear();
      
      if(m_totalRateShapingBuf == 0 && m_complete == true && m_srcSsrcSet.empty())
      {
        NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"All  streams are completed");
        Simulator::Cancel(m_sendEvent);
        Simulator::Cancel(m_enqueueEvent);
        
        if(!m_receiving)
        {
          NS_LOG_INFO("RTCP Stop when send BYE "<<m_localSsrc);
          Simulator::Cancel(m_rtcpEvent);
          Simulator::Cancel(m_rembEvent);
          SendRtcp(*header, false);
          return;
        }
      }
    }

    SendRtcp(*header, true);
}

void GccNode::SendRtcp(GccRtcpHeader header, bool reschedule)
{
    auto packet = Create<Packet> ();
    packet->AddHeader (header);
    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::SendRtcp, " << packet->ToString ());
    m_socket->SendTo (packet, 0, InetSocketAddress{m_destIp, m_destPort});

    m_sending = false;
    m_receiving = false;
   
    if(reschedule)
    { 
      m_rtcpEvent = Simulator::Schedule(ns3::MilliSeconds(GetNextRtcpTime()), &GccNode::CreateRtcp, this);
    }
}

void GccNode::SendRemb()
{
  RembHeader header(RembHeader::RTP_REMB);
  NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::SendRemb");

  header.SetSendSsrc(m_localSsrc);

  std::set<uint32_t> rembSet = m_recvSsrcSet;
  rembSet.insert(m_localSsrc);
  NS_ASSERT(rembSet.size() <= RembHeader::MAX_SSRCS_NUM);

  uint32_t bitrate = m_recvController->GetBitrate();; 
  NS_ASSERT(header.AddRembFeedback(bitrate, rembSet) == RembHeader::RTCP_NONE);

  Ptr<Packet> packet = Create<Packet> ();
  packet->AddHeader (header);

  NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::SendRemb, " << packet->ToString ());
  m_socket->SendTo (packet, 0, InetSocketAddress{m_destIp, m_destPort});

  m_rembEvent = Simulator::Schedule (ns3::MilliSeconds(REMBINTERVAL), &GccNode::SendRemb, this); //REMB is sent per 1s (webRTC)
}

uint32_t GccNode::GetNextRtcpTime()
{
    uint32_t minIntervalMs = 0;  
    // Calculate bandwidth for video; 360s / send bandwidth in kbit/s. (WebRTC)
    uint32_t send_bitrate_kbit = m_rSend / 1000;
         
    if (send_bitrate_kbit != 0)
      minIntervalMs = 360*1000 / send_bitrate_kbit;

    if (minIntervalMs > VIDEOINTERVAL)  //assume that we send video... if audio, use AUDIOINTERVAL
    {     
      minIntervalMs = VIDEOINTERVAL;
    }
      
    // The interval between RTCP packets is varied randomly over the
    // range [1/2,3/2] times the calculated interval.
    uint32_t rtcpNext = rand()%(minIntervalMs+1) + minIntervalMs/2;

    NS_LOG_INFO("Next RTCP Send Interval : "<<rtcpNext);

    return rtcpNext;
}

/*Receive Methods*/
void GccNode::RecvPacket (Ptr<Socket> socket)
{

    Address remoteAddr;
    auto Packet = m_socket->RecvFrom (remoteAddr);
    NS_ASSERT (Packet);
    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::RecvPacket, " << Packet->ToString ());

    auto rIPAddress = InetSocketAddress::ConvertFrom (remoteAddr).GetIpv4 ();
    auto rport = InetSocketAddress::ConvertFrom (remoteAddr).GetPort ();
    NS_ASSERT (rIPAddress == m_destIp);
    NS_ASSERT (rport == m_destPort);

    //check payload type
    auto tmp_packet = Packet->Copy();
    RtpHeader tmp_header{};
    tmp_packet->RemoveHeader(tmp_header);
    uint32_t pt;
    if(tmp_header.IsMarker()) //if RTCP packet, marker is true otherwise 0
    {
      auto tmp_packet2 = Packet->Copy();
      RtcpHeader tmp_header2{};
      tmp_packet2->RemoveHeader(tmp_header2);
      pt = tmp_header2.GetPacketType();
    }
    else
      pt = tmp_header.GetPayloadType();

    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::packet type : "<<pt);
    switch(pt)
    {
      case 96:
        RecvDataPacket(Packet, remoteAddr);
        break;
      case GccRtcpHeader::RTCP_SR:
      case GccRtcpHeader::RTCP_RR:
        RecvRtcp(Packet, remoteAddr);
        break;
      case GccRtcpHeader::RTP_REMB:
        RecvRembPacket(Packet, remoteAddr);
        break;
      default:
        NS_LOG_DEBUG("Invalid Packet is received");
        exit(1);
    }

   // m_rSend = m_senderController->getBitrate();
}

void GccNode::RecvDataPacket(Ptr<Packet> p, Address remoteAddr)
{

    uint64_t nowMs = Simulator::Now().GetMilliSeconds();
    uint32_t pSize = p->GetSize();

    RtpHeader header{};
    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::RecvDataPacket, " << p->ToString ());
    p->RemoveHeader (header);

    uint32_t recvSsrc = header.GetSsrc();
    uint32_t seq = header.GetSequence();

    if(m_recvSsrcSet.find(recvSsrc) == m_recvSsrcSet.end())
    {
      m_recvSsrcSet.insert(recvSsrc);
      m_recvSeq[recvSsrc] = seq-1;
    }
    
    m_lost[recvSsrc] += seq - (m_recvSeq[recvSsrc]+1);
    m_cumLost[recvSsrc] += m_lost[recvSsrc];
    m_recvPackets[recvSsrc]++;
    m_recvSeq[recvSsrc] = seq;

    m_recvController->UpdateDelayBasedBitrate(nowMs, seq, header.GetTimestamp(), nowMs, p->GetSize(), 0);

    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::RecvDataPacket, Lost : "<<m_lost[recvSsrc]<<", cumLost : "<<m_cumLost[recvSsrc]<<", recvPackets : "<<m_recvPackets[recvSsrc]++);

    m_receiving = true;

    //Print Log
    uint32_t delay = uint32_t(nowMs)-header.GetTimestamp();
    m_pDelay = m_pDelay*0.7 + delay*0.3; //ms
    m_recvDataBytes += pSize;
    m_recvAllBytes += pSize;

    if(ns3::MilliSeconds(nowMs) > m_lastThroCheckTime + ns3::MilliSeconds(THROCHECKINTERVAL))
    {
      NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<"\tGoodput : "<<m_recvDataBytes*8*1000/((double)THROCHECKINTERVAL*1000*1000)<<"\tThroughput : "<<m_recvAllBytes*8*1000/((double)THROCHECKINTERVAL*1000*1000));
      m_lastThroCheckTime = m_lastThroCheckTime + ns3::MilliSeconds(THROCHECKINTERVAL);
      m_recvDataBytes = 0;
      m_recvAllBytes = 0;
    }

    if(ns3::MilliSeconds(nowMs) > m_lastDelayCheckTime + ns3::MilliSeconds(DELAYCHECKINTERVAL))
    {
      NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<"\tavg pDelay : "<<m_pDelay<<"\tinstant pDelay : "<<delay);
      m_lastDelayCheckTime = m_lastDelayCheckTime + ns3::MilliSeconds(DELAYCHECKINTERVAL);
    }
    
}

void GccNode::RecvRtcp(Ptr<Packet> p, Address remoteAddr)
{
    uint64_t nowMs = Simulator::Now().GetMilliSeconds();
    m_recvAllBytes += p->GetSize();
    GccRtcpHeader header{};
    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::RecvRtcp, " << p->ToString ());
    p->RemoveHeader (header);

    if(m_remoteSsrc == 0)
      m_remoteSsrc = header.GetSendSsrc();
    else
      NS_ASSERT(m_remoteSsrc == header.GetSendSsrc());

    GccRtcpHeader::SendReportBlock srb;
    if(header.GetPacketType() == GccRtcpHeader::RTCP_SR)
      GccRtcpHeader::SendReportBlock srb = header.GetSRB();

    //add action for SR Block, if required...
    std::vector<GccRtcpHeader::RecvReportBlock> rrbs{};
    rrbs = header.GetRRBs ();

    GccRtcpHeader::ByeBlock bye = header.GetBye();

    if(!rrbs.empty())
    {
      for(const auto& rrb : rrbs)
      {
        NS_ASSERT(m_localSsrc == rrb.m_sourceSsrc || (m_srcSsrcSet.find(rrb.m_sourceSsrc) != m_srcSsrcSet.end()));
     
        uint32_t rtt = nowMs-rrb.m_lastSRTime-rrb.m_SRDelay;
     
        NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::RecvRrBlocks, ssrc : "<<rrb.m_sourceSsrc<<", frac : "<<rrb.m_fractionLost<<", cumLost : "<<rrb.m_cumNumLost<<", seq : "<<rrb.m_highestSeqNum<<", last SR Recv Time : "<<rrb.m_lastSRTime<<", SR Delay : "<<rrb.m_SRDelay<<", rtt : "<<rtt);
        
        m_lastSrRecvTime[rrb.m_sourceSsrc] = Simulator::Now().GetMilliSeconds();
      }
      m_senderController->OnReceivedRtcpReceiverReportBlocks(rrbs, nowMs);
    }
    
    if(!bye.m_ssrcSet.empty())
    {
      NS_ASSERT(bye.m_ssrcSet.size() > 0);
      
      for(auto const& b_ssrc : bye.m_ssrcSet)
      {
        NS_ASSERT(m_recvSsrcSet.find(b_ssrc) != m_recvSsrcSet.end());
        m_recvSsrcSet.erase(b_ssrc);
        m_recvSeq.erase(b_ssrc);
        m_lost.erase(b_ssrc);
        m_cumLost.erase(b_ssrc);
        m_recvPackets.erase(b_ssrc);

        if(m_recvSsrcSet.empty())
        {
          Simulator::Cancel (m_rembEvent);
          if(m_complete==true && m_srcSsrcSet.empty())
          {
            NS_LOG_INFO("RTCP Stop "<<m_localSsrc);
            Simulator::Cancel(m_rtcpEvent);
          }
        }
      }
    }
    return;
}

void GccNode::RecvRembPacket(Ptr<Packet> p, Address remoteAddr)
{
    m_recvAllBytes += p->GetSize();
    RembHeader header{};
    NS_LOG_INFO(Simulator::Now().ToDouble(Time::S)<<" "<<"GccNode::RecvRembPacket, " << p->ToString ());
    p->RemoveHeader (header);

    RembHeader::RembBlock remb = header.GetRembBlock();
    //gcc feedback..
}


}
