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
#define INITRATE 1000000 //1Mbps
#define VIDEOINTERVAL 1000000 //1s == 1000000us
#define AUDIOINTERVAL 5000000 //1s == 1000000us
#define REMBINTERVAL 1000000 //1s == 1000000us

namespace ns3 {

GccNode::GccNode ()
: m_destIp{}
, m_destPort{0}
, m_localPort{0}
, m_socket{NULL}
, m_localSsrc{rand()}
, m_numSsrc{1}
, m_enqueueEvent{}
, m_sendEvent{}
, m_sendOversleepEvent{}
, m_rtcpEvent{}
, m_rembEvent{}
, m_rSend{INITRATE}
, m_totalRateShapingBuf{0}
, m_nextSendTstmpUs{0}
, m_nextRtcpTstmpUs{0}
, m_nextEnqSsrcIndex{0}
, m_nextSendSsrcIndex{0}
, m_sending{false}
, m_remoteSsrc{0}
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
//    m_sendController = std::make_shared<rmcat::GccSendController> ();

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

  m_numSsrc += m_srcSsrcSet.size();

  return true;
}

void GccNode::StartApplication ()
{
    // RTP initial values for sequence number and timestamp SHOULD be random (RFC 3550)
    m_sequence[m_localSsrc] = rand ();
    m_rtpTsOffset[m_localSsrc] = rand ();
    m_rateShapingBytes[m_localSsrc] = 0;
    m_enqBytes[m_localSsrc] = 0;
    m_complete[m_localSsrc] = false;


    for(const auto& ssrc : m_srcSsrcSet)
    {
      m_sequence[ssrc] = rand();
      m_rtpTsOffset[ssrc] = rand();
      m_rateShapingBytes[ssrc] = 0;
      m_enqBytes[ssrc] = 0;
      m_complete[ssrc] = false;
    }

    m_rSend = INITRATE;

    m_nextEnqSsrcIndex = 0;
    m_nextSendSsrcIndex = 0;
    m_nextSendTstmpUs = 0;

    m_enqueueEvent = Simulator::Schedule (Seconds (0.0), &GccNode::EnqueuePacket, this);
    
    uint32_t minIntervalUs = 0;  
    // Calculate bandwidth for video; 360s / send bandwidth in kbit/s. (WebRTC)
    uint32_t send_bitrate_kbit = m_rSend / 1000;
         
    if (send_bitrate_kbit != 0)
      minIntervalUs = 360000 / send_bitrate_kbit;

    if (minIntervalUs > VIDEOINTERVAL)  //assume that we send video... if audio, use AUDIOINTERVAL
    {     
      minIntervalUs = VIDEOINTERVAL;
    }
      
    // The interval between RTCP packets is varied randomly over the
    // range [1/2,3/2] times the calculated interval.
    uint32_t rtcpNext = rand()%(minIntervalUs+1) + minIntervalUs/2;
    const uint64_t nowUs = Simulator::Now ().GetMicroSeconds ();

    m_nextRtcpTstmpUs = nowUs + rtcpNext;

    m_rtcpEvent = Simulator::Schedule(ns3::MicroSeconds(m_nextRtcpTstmpUs), &GccNode::SendSr, this);
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

    for(uint32_t i = 0;i<m_numSsrc;i++)
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

      if(m_maxSize[enq_ssrc] != 0 && m_enqBytes[enq_ssrc] >= m_maxSize[enq_ssrc])
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

    NS_LOG_INFO ("GccNode::EnqueuePacket, packet enqueued, packet length: " << bytesToSend
                 << ", buffer size: " << m_rateShapingBuf[enq_ssrc].size ()
                 << ", buffer bytes: " << m_rateShapingBytes[enq_ssrc] << ", total enqueue bytes : "<<m_enqBytes[enq_ssrc]);

    double secsToNextEnqPacket = codec->second;
    Time tNext{Seconds (secsToNextEnqPacket)};
    m_enqueueEvent = Simulator::Schedule (tNext, &GccNode::EnqueuePacket, this);

    if (m_totalRateShapingBuf == 1) {
        // Buffer was empty
        const uint64_t nowUs = Simulator::Now ().GetMicroSeconds ();
        const uint64_t usToNextSentPacket = nowUs < m_nextSendTstmpUs ?
                                                    m_nextSendTstmpUs - nowUs : 0;
        NS_LOG_INFO ("(Re-)starting the send timer: nowUs " << nowUs
                     << ", bytesToSend " << bytesToSend
                     << ", usToNextSentPacket " << usToNextSentPacket
                     << ", m_rSend " << m_rSend
                     << ", secsToNextEnqPacket " << secsToNextEnqPacket);

        Time tNext{MicroSeconds (usToNextSentPacket)};
        m_sendEvent = Simulator::Schedule (tNext, &GccNode::SendPacket, this, usToNextSentPacket);
    }
}

void GccNode::SendPacket (uint64_t usSlept)
{
    uint32_t send_ssrc = 0;
    bool can_send = false;
    for(uint32_t i=0;i<m_numSsrc;i++)
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
        bool allComplete = true;
        if(m_maxSize[send_ssrc] != 0 && (m_enqBytes[send_ssrc] >= m_maxSize[send_ssrc] && !m_complete[send_ssrc]))
        {
          m_complete[send_ssrc] = true;
          m_byeSet.insert(send_ssrc);
          m_srcSsrcSet.erase(send_ssrc);
          m_numSsrc--;
          for(auto const& c : m_complete)
          {
            if(!c.second)
            {
              allComplete = false;
              break;
            }
          }
        }
        
        if(allComplete)
        {
          NS_LOG_INFO("All stream is completed");
          m_sending = false;
          StopApplication();
          return;
        }

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
      NS_LOG_INFO("All stream buffer is empty but not completed.");
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

    NS_LOG_INFO ("GccNode::SendPacket, packet dequeued, packet length: " << bytesToSend
                 << ", buffer size: " << m_rateShapingBuf[send_ssrc].size ()
                 << ", buffer bytes: " << m_rateShapingBytes[send_ssrc]);

    // Synthetic oversleep: random uniform [0% .. 1%]
    uint64_t oversleepUs = usSlept * (rand () % 100) / 10000;
    Time tOver{MicroSeconds (oversleepUs)};
    m_sendOversleepEvent = Simulator::Schedule (tOver, &GccNode::SendOverSleep,
                                                this, bytesToSend,send_ssrc);

    // schedule next sendData
    const double usToNextSentPacketD = double (bytesToSend) * 8. * 1000. * 1000. / m_rSend;
    const uint64_t usToNextSentPacket = uint64_t (usToNextSentPacketD);

    if (m_totalRateShapingBuf == 0) {
        // Buffer became empty
        const auto nowUs = Simulator::Now ().GetMicroSeconds ();
        m_nextSendTstmpUs = nowUs + usToNextSentPacket;
        return;
    }

    Time tNext{MicroSeconds (usToNextSentPacket)};
    m_sendEvent = Simulator::Schedule (tNext, &GccNode::SendPacket, this, usToNextSentPacket);
}

void GccNode::SendOverSleep (uint32_t bytesToSend, uint32_t send_ssrc) 
{
    if(!m_sending)
      m_sending = true;
    
    const auto nowUs = Simulator::Now ().GetMicroSeconds ();
    //m_sendController->... required

    ns3::RtpHeader header{96}; // 96: dynamic payload type, according to RFC 3551
    header.SetSequence (m_sequence[send_ssrc]++);
    NS_ASSERT (nowUs >= 0);
    // Most video payload types in RFC 3551, Table 5, use a 90 KHz clock
    // Therefore, assuming 90 KHz clock for RTP timestamps
    header.SetTimestamp (m_rtpTsOffset[send_ssrc] + uint32_t (nowUs * 90 / 1000));;
    header.SetSsrc (send_ssrc);

    auto packet = Create<Packet> (bytesToSend);
    packet->AddHeader (header);

    NS_LOG_INFO ("GccNode::SendOverSleep, " << packet->ToString ());
    m_socket->SendTo (packet, 0, InetSocketAddress{m_destIp, m_destPort});
}


void GccNode::SendSr()
{
    if(!m_sending)
    {
      SendRr();
      return;
    }

    GccRtcpHeader* header = new GccRtcpHeader(GccRtcpHeader::RTCP_SR);
    header->SetSendSsrc(m_localSsrc);

    NS_ASSERT(header->AddSRFeedback(0,0,0,0,0) == GccRtcpHeader::RTCP_NONE); //Data in SR does not affect Gcc Performance

    uint32_t rrSize = m_srcSsrcSet.size();

    for(const auto& ssrc : m_recvSsrcSet)
    {
      double fractionRatio = m_lost[ssrc]/(double)m_recvPackets[ssrc];
      uint8_t frac = uint8_t(255.0 * fractionRatio);
      
      uint32_t ret = 0;
      ret = header->AddRRFeedback(ssrc, frac, m_cumLost[ssrc], m_recvSeq[ssrc], 0, 0, 0);
      if(ret == GccRtcpHeader::RTCP_TOO_LONG)
      {
        header->SetTypeOrCount(31); //maximum number of RR blocks by RFC & WebRTC
        rrSize -= 31;
        SendRtcp(*header, false);
        
        delete header;
        header = new GccRtcpHeader(GccRtcpHeader::RTCP_RR);
      }
      NS_ASSERT(ret == GccRtcpHeader::RTCP_NONE);
    }

    header->SetTypeOrCount(rrSize);

    if(!m_byeSet.empty())
    { 
      NS_ASSERT(header->AddByeFeedback(m_byeSet) == GccRtcpHeader::RTCP_NONE);
      m_byeSet.clear();
    }

    SendRtcp(*header, true);
}

void GccNode::SendRr()
{
    GccRtcpHeader* header = new GccRtcpHeader(GccRtcpHeader::RTCP_RR);
   
    header->SetSendSsrc(m_localSsrc);
    
    uint32_t rrSize = m_srcSsrcSet.size();
    
    for(const auto& ssrc : m_recvSsrcSet)
    {
      double fractionRatio = m_lost[ssrc]/(double)m_recvPackets[ssrc];
      uint8_t frac = uint8_t(255.0 * fractionRatio);
      uint32_t ret = 0;

      ret = header->AddRRFeedback(ssrc, frac, m_cumLost[ssrc], m_recvSeq[ssrc], 0, 0, 0);
      if(ret == GccRtcpHeader::RTCP_TOO_LONG)
      {
        header->SetTypeOrCount(31); //maximum number of RR blocks by RFC & WebRTC
        rrSize -= 31;
        SendRtcp(*header, false);
        
        delete header;
        header = new GccRtcpHeader(GccRtcpHeader::RTCP_RR);
      }
      
      NS_ASSERT(ret == GccRtcpHeader::RTCP_NONE);
      m_lost[ssrc] = 0;
      m_recvPackets[ssrc] = 0;
    }
    
    header->SetTypeOrCount(rrSize);

    if(!m_byeSet.empty()) 
    {
      NS_ASSERT(header->AddByeFeedback(m_byeSet) == GccRtcpHeader::RTCP_NONE);
      m_byeSet.clear();
    }

    SendRtcp(*header, true);
}

void GccNode::SendRtcp(GccRtcpHeader header, bool reschedule)
{
    auto packet = Create<Packet> ();
    packet->AddHeader (header);
    NS_LOG_INFO ("GccNode::SendRtcp, " << packet->ToString ());
    m_socket->SendTo (packet, 0, InetSocketAddress{m_destIp, m_destPort});

    if(reschedule)
    {
      Time tNext {MicroSeconds (0)}; //should set!
      m_rtcpEvent = Simulator::Schedule (tNext, &GccNode::SendSr, this);
    }

    m_sending = false;
}

void GccNode::SendRemb()
{
  RembHeader header(RembHeader::RTP_REMB);

  header.SetSendSsrc(m_localSsrc);

  std::set<uint32_t> rembSet = m_recvSsrcSet;
  rembSet.insert(m_localSsrc);
  NS_ASSERT(rembSet.size() <= RembHeader::MAX_SSRCS_NUM);

  uint32_t bitrate = 0; //set..
  NS_ASSERT(header.AddRembFeedback(bitrate, rembSet) == RembHeader::RTCP_NONE);

  auto packet = Create<Packet> ();
  packet->AddHeader (header);
  NS_LOG_INFO ("GccNode::SendRemb, " << packet->ToString ());
  m_socket->SendTo (packet, 0, InetSocketAddress{m_destIp, m_destPort});

  m_rembEvent = Simulator::Schedule (ns3::MicroSeconds(REMBINTERVAL), &GccNode::SendRemb, this); //REMB is sent per 1s (webRTC)
}

/*Receive Methods*/
void GccNode::RecvPacket (Ptr<Socket> socket)
{
    if(!m_rembEvent.IsRunning())
      m_rembEvent = Simulator::Schedule (ns3::MicroSeconds(REMBINTERVAL), &GccNode::SendRemb, this); //REMB is sent per 1s (webRTC)

    Address remoteAddr;
    auto Packet = m_socket->RecvFrom (remoteAddr);
    NS_ASSERT (Packet);

    auto rIPAddress = InetSocketAddress::ConvertFrom (remoteAddr).GetIpv4 ();
    auto rport = InetSocketAddress::ConvertFrom (remoteAddr).GetPort ();
    NS_ASSERT (rIPAddress == m_destIp);
    NS_ASSERT (rport == m_destPort);

    //check payload type
    auto tmp_packet = Packet->Copy();
    RtpHeader tmp_header{};
    tmp_packet->RemoveHeader(tmp_header);

    uint32_t pt = tmp_header.GetPayloadType();
    switch(pt)
    {
      case 96:
        RecvDataPacket(Packet, remoteAddr);
        break;
      case GccRtcpHeader::RTCP_SR:
        RecvSrPacket(Packet, remoteAddr);
        break;
      case GccRtcpHeader::RTCP_RR:
        RecvRrPacket(Packet, remoteAddr);
        break;
      case GccRtcpHeader::RTP_REMB:
        RecvRembPacket(Packet, remoteAddr);
        break;
      default:
        NS_LOG_DEBUG("Invalid Packet is received");
        exit(1);
    }
}

void GccNode::RecvDataPacket(Ptr<Packet> p, Address remoteAddr)
{
    RtpHeader header{};
    NS_LOG_INFO ("GccNode::RecvDataPacket, " << p->ToString ());
    p->RemoveHeader (header);

    uint32_t recvSsrc = header.GetSsrc();
    uint32_t seq = header.GetSequence();

    if(m_recvSsrcSet.find(recvSsrc) == m_recvSsrcSet.end())
      m_recvSsrcSet.insert(recvSsrc);
    
    m_lost[recvSsrc] += seq - m_recvSeq[recvSsrc];
    m_cumLost[recvSsrc] += m_lost[recvSsrc];
    m_recvPackets[recvSsrc]++;

    //AddFeedback (header.GetSequence (), recvTimestampUs);
}

void GccNode::RecvSrPacket(Ptr<Packet> p, Address remoteAddr)
{
    GccRtcpHeader header{};
    NS_LOG_INFO ("GccNode::RecvSrPacket, " << p->ToString ());
    p->RemoveHeader (header);

    if(m_remoteSsrc == 0)
    {
      NS_ASSERT(m_recvSsrcSet.find(header.GetSendSsrc()) != m_recvSsrcSet.end());
      m_remoteSsrc = header.GetSendSsrc();
    } 
    else
      NS_ASSERT(m_remoteSsrc == header.GetSendSsrc());

    GccRtcpHeader::SendReportBlock srb = header.GetSRB();
    //add action for SR Block, if required...

    std::vector<GccRtcpHeader::RecvReportBlock> rrbs{};
    rrbs = header.GetRRBs ();

    GccRtcpHeader::ByeBlock bye = header.GetBye();

    if(!rrbs.empty())
    {
      for(const auto& rrb : rrbs)
      {
        NS_ASSERT(m_srcSsrcSet.find(rrb.m_sourceSsrc) != m_srcSsrcSet.end());
      }
      //feedback to controller
    }
    
    if(!bye.m_ssrcSet.empty())
    {
      NS_ASSERT(bye.m_ssrcSet.size() > 0);
      
      for(auto const& b : bye.m_ssrcSet)
      {
        NS_ASSERT(m_recvSsrcSet.find(b) != m_recvSsrcSet.end());
        m_recvSsrcSet.erase(b);
      }
    }
    return;
}

void GccNode::RecvRrPacket(Ptr<Packet> p, Address remoteAddr)
{
    NS_LOG_INFO ("GccNode::RecvRrPacket, " << p->ToString ());
    
    GccRtcpHeader header{};
    p->RemoveHeader (header);
    
    if(m_remoteSsrc == 0)
    {
      NS_ASSERT(m_recvSsrcSet.find(header.GetSendSsrc()) != m_recvSsrcSet.end());
      m_remoteSsrc = header.GetSendSsrc();
    } 
    else
      NS_ASSERT(m_remoteSsrc == header.GetSendSsrc());
    
    std::vector<GccRtcpHeader::RecvReportBlock> rrbs{};
    rrbs = header.GetRRBs ();
    GccRtcpHeader::ByeBlock bye = header.GetBye();

    if(!rrbs.empty())
    {
      for(const auto& rrb : rrbs)
      {
        NS_ASSERT(m_srcSsrcSet.find(rrb.m_sourceSsrc) != m_srcSsrcSet.end());
      }
      //feedback to controller
    }
    
    if(!bye.m_ssrcSet.empty())
    {
      NS_ASSERT(bye.m_ssrcSet.size() > 0);
      
      for(auto const& b : bye.m_ssrcSet)
      {
        NS_ASSERT(m_recvSsrcSet.find(b) != m_recvSsrcSet.end());
        m_recvSsrcSet.erase(b);
      }
    }
    
    return;
}

void GccNode::RecvRembPacket(Ptr<Packet> p, Address remoteAddr)
{
    RembHeader header{};
    NS_LOG_INFO ("GccNode::RecvSrPacket, " << p->ToString ());
    p->RemoveHeader (header);

    RembHeader::RembBlock remb = header.GetRembBlock();
    //gcc feedback..
}


}
