/******************************************************************************
 * q
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
 * Gcc Sender Side controller (Loss Based Congestion Control) implementation for gcc ns3module.
 *
 * @version 0.1.0
 */
#include "gcc-sender-controller.h"
#include <sstream>
#include <cassert>

namespace rmcat {
const int kDefaultMaxBitrateBps = 1000000000;
const int64_t kLowBitrateLogPeriodMs = 10000;
const int kMinBitrateBps = 10000;
const int kMaxBitrateBps = 1500000;
const int kInitialBitrateBps = 300000;

GccSenderController::GccSenderController() :

    m_lastTimeCalcUs{0},
    m_lastTimeCalcValid{false},
    m_QdelayUs{0},
    m_ploss{0},
    m_plr{0.f},
    min_configured_bitrate_bps_(150000),
    max_configured_bitrate_bps_(1500000),
    
    m_sending_rate_(max_configured_bitrate_bps_),
    
    min_bitrate_configured_{min_configured_bitrate_bps_},
    max_bitrate_configured_{max_configured_bitrate_bps_},
    last_low_bitrate_log_ms_(-1),
    m_Bitrate_valid_(false){}

GccSenderController::~GccSenderController() {}

/*
void GccSenderController::setCurrentBw(float newBw) {
    m_initBw = newBw;
}
*/

void GccSenderController::reset() {
    m_lastTimeCalcUs = 0;
    m_lastTimeCalcValid = false;

    m_QdelayUs = 0;
    m_ploss = 0;
    m_plr = 0.f;
}

float getBitrate() {
    return m_sending_rate_;
}

void GccSenderController::SetSendBitrate(int bitrate, int nowms) {
    if (bitrate > max_bitrate_configured_){
         bitrate = max_bitrate_configured_;
    }

    if (bitrate < min_bitrate_configured_){
        if(last_low_bitrate_log_ms_ == -1 || now_ms - last_low_bitrate_log_ms > kLowBitrateLogPeriodMs) {
	    last_low_bitrate_log_ms_ = now_ms;
	}

	bitrate = min_bitrate_configured_;
    }    

    m_sending_rate_ = bitrate;
}

void GccSenderController::SetMinMaxBitrate(int min_bitrate, int max_bitrate) {
    min_bitrate_configured_ = std::max(min_bitrate, min_bitrate_configured_);

    if(max_bitrate > 0) {
	max_bitrate_configured_ = std::max<uint32_t>(min_bitrate_configured_, max_bitrate);
    } else {
	max_bitrate_configured_ = kDefaultMaxBitrateBps;
    }
}

void GccSenderController::ApplyDelayBasedBitrate(float DelayBasedEstimateBitrate) {
    // TODO Called when the REMB messages are received from receiver.
    // Calculate proper m_sending_rate_ according to delay-based estimated bitrate.
    
    if (!m_Bitrate_valid_){
        // First ApplyDelayBasedBitrate called. We need to initialize m_sending_rate_ and min, max configured bitrate.
        m_Bitrate_valid_ = true;

	SetMinMaxBitrate(kMinBitrateBps, kMaxBitrateBps);
	SetSendBitrate(ns3::Simulator::now().GetMilliSeconds(), kInitialBitrateBps);

	return;
    }
    if (DelayBasedEstimateBitrate > m_sending_rate_)        DelayBasedEstimateBitrate = m_sending_rate_;

    if (DelayBasedEstimateBitrate > max_bitrate_configured_) {
        DelayBasedEstimateBitrate = max_bitrate_configured_;
    }

    if (DelayBasedEstimateBitrate < min_bitrate_configured_) {
        if (last_low_bitrate_log_ms_ == -1 || now_ms - last_low_bitrate_log_ms_ > kLowBitrateLogPeriodMs) {
            last_low_bitrate_log_ms_ = now_ms;
	}

    	DelayBasedEstimateBitrate = min_bitrate_configured_;
    }

    m_sending_rate_ = DelayBasedEstimateBitrate;
}


void GccSenderController::UpdateLossBasedBitrate(uint64_t nowMs,
                                                 float plr) {
    // TODO Loss Based Controller Implementation
   
}

/*
float GccSenderController::getBandwidth(uint64_t nowUs) const {

    return m_initBw;
}
*/

/*
void GccSenderController::updateMetrics() {
    uint64_t qdelayUs;
    bool qdelayOK = getCurrentQdelay(qdelayUs);
    if (qdelayOK) m_QdelayUs = qdelayUs;

    float rrate;
    bool rrateOK = getCurrentRecvRate(rrate);
    if (rrateOK) m_RecvR = rrate;

    uint32_t nLoss;
    float plr;
    bool plrOK = getPktLossInfo(nLoss, plr);
    if (plrOK) {
        m_ploss = nLoss;
        m_plr = plr;
    }
}
*/

void GccSenderController::logStats(uint64_t nowUs) const {

    std::ostringstream os;
    os << std::fixed;
    os.precision(RMCAT_LOG_PRINT_PRECISION);

    os  << " algo:dummy " << m_id
        << " ts: "     << (nowUs / 1000)
        << " loglen: " << m_packetHistory.size()
        << " qdel: "   << (m_QdelayUs / 1000)
        << " ploss: "  << m_ploss
        << " plr: "    << m_plr
        << " rrate: "  << m_RecvR
        << " srate: "  << m_initBw;
    logMessage(os.str());
}

}
