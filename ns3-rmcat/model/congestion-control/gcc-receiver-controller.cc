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
 * Gcc Receiver Side Controller (Loss Based Congestion Control) implementation for gcc ns3 module.
 *
 * @version 0.1.0
 */
#include "gcc-receiver-controller.h"
#include <sstream>
#include <cassert>
#define BURST_TIME 5

namespace rmcat {

GccRecvController::GccRecvController() :
    m_lastTimeCalcUs{0},
    m_lastTimeCalcValid{false},
    m_QdelayUs{0},
    m_ploss{0},
    m_plr{0.f},
    m_RecvR{0.}

    curr_group_num_{0}
    curr_group_sseq_{0}
    curr_group_stime_{0}
    curr_group_eseq_{0}
    curr_group_etime_arrival_{0}
    curr_group_etime_departure_{0}

    prev_group_sseq_{0}
    prev_group_stime_{0}
    prev_group_eseq_{0}
    prev_group_etime_arrival_{0}
    prev_group_etime_departure_{0}

    curr_group_size_{0}
    prev_group_size_{0}
    group_size_interval_{0}

    prev_pkt_seq_{0}
    prev_pkt_txTime_{0}
    prev_pkt_rxTime_{0}

    i_arrival_{0}
    i_departure_{0}
    i_delay_var_{0}

    m_group_changed_{false}
    m_first_packet_{false}{}
    

GccRecvController::~GccRecvController() {}

void GccRecvController::setCurrentBw(float newBw) {
    m_initBw = newBw;
}

void GccRecvController::reset() {
    m_lastTimeCalcUs = 0;
    m_lastTimeCalcValid = false;

    m_QdelayUs = 0;
    m_ploss = 0;
    m_plr = 0.f;
    m_RecvR = 0.;

    curr_group_num_ = 0;
}

// txTimestamp : current received packet's txTimestamp in rtp header.
// rxTimestamp : current received packet's rxTimestamp. (packet received time in local)
// sequence : current received packet's sequence.
void GccRecvController::UpdateGroupInfo(uint64_t nowUs, uint16_t sequence, uint64_t txTimestamp,
					uint64_t rxTimestamp, uint64_t packet_size){
    m_group_changed_ = false;

    if(m_first_packet_) {
	// First Packet... We need to initialize variable.
        curr_group_num_ = 0;
	curr_group_sseq_ = sequence;
	curr_group_stime_ = txTimestamp;
	
	curr_group_size_ = 0;
        curr_group_size_ += packet_size;
	
	prev_pkt_seq_ = sequence;
	prev_pkt_txTime_ = txTimestamp;
	prev_pkt_rxTime_ = rxTimestamp;
        
	m_first_packet_ = false;
	return;
    }

    if (curr_group_stime_ + BURST_TIME < txTimestamp) {
        // Groups are change.
	if(curr_group_num_ == 0) {
            // First group ends... There's no inter variables calculation.
	    // Save current group information to prev_group_*.
	    prev_group_sseq_ = curr_group_sseq_;
	    prev_group_stime_ = curr_group_stime_;
	    prev_group_eseq_ = prev_pkt_seq_;
	    prev_group_etime_arrival_ = prev_pkt_rxTime_;
	    prev_group_etime_departure_ = prev_pkt_txTime_;
            
	    prev_group_size_ = curr_group_size_;
	    curr_group_size_ = packet_size;
	    
	    // Update new group information.
	    curr_group_num_++;
	    curr_group_sseq_ = sequence;
	    curr_group_stime_ = txTimestamp;
	}
	else {
	    /* curr_group_etime_arrival_ = prev_pkt_seq_rxTime_;
	    curr_group_etime_departure_ = prev_pkt_seq_txTime_; */
	    i_arrival_ = prev_pkt_rxTime_ - prev_group_etime_arrival_;
	    i_departure_ = prev_pkt_txTime_ - prev_group_etime_departure_;

	    i_delay_var_ = i_arrival_ - i_departure_;

	    // Save Current group information into prev_group_*.
            prev_group_sseq_ = curr_group_sseq_;
	    prev_group_stime_ = curr_group_stime_;
	    prev_group_eseq_ = prev_pkt_seq_;
	    prev_group_etime_arrival_ = prev_pkt_rxTime_;
	    prev_group_etime_departure_ = prev_pkt_txTime_;

	    group_size_interval_ = prev_group_size_ - curr_group_size_;
	    prev_group_size_ = curr_group_size_;
	    curr_group_size_ = 0;

	    //Update new group information.
	    curr_group_num++;
	    curr_group_sseq_ = sequence;
	    curr_group_stime_ = txTimestamp;

	    m_group_changed_ = true;
	}
        
	curr_group_size_ += packet_size;

	prev_pkt_seq_ = sequence;
	prev_pkt_txTime_ = txTimestamp;
	prev_pkt_rxTime_ = rxTimestamp:;

	return;
    }

    int temp_arrival = prev_pkt_rxTime_ - prev_group_etime_arrival_;
    int temp_departure = prev_pkt_txTime_ - prev_group_etime_departure_;
    int temp_delay_var = temp_arrival - temp_departure;

    if(temp_arrival < BURST_TIME && temp_delay_var < 0) {
        curr_group_size_ += packet_size;
	
	prev_pkt_seq_ = sequence;
	prev_pkt_txTime_ = txTimestamp;
	prev_pkt_rxTime_ = rxTimestamp;

	return;
    }
    else {
        // Groups are change.
	if(curr_group_num_ == 0) {
            // First group ends... There's no inter variables calculation.
	    // Save current group information to prev_group_*.
	    prev_group_sseq_ = curr_group_sseq_;
	    prev_group_stime_ = curr_group_stime_;
	    prev_group_eseq_ = prev_pkt_seq_;
	    prev_group_etime_arrival_ = prev_pkt_rxTime_;
	    prev_group_etime_departure_ = prev_pkt_txTime_;
            
	    prev_group_size_ = curr_group_size_;
	    curr_group_size_ = packet_size;
	    
	    // Update new group information.
	    curr_group_num_++;
	    curr_group_sseq_ = sequence;
	    curr_group_stime_ = txTimestamp;
	}
	else {
	    /* curr_group_etime_arrival_ = prev_pkt_seq_rxTime_;
	    curr_group_etime_departure_ = prev_pkt_seq_txTime_; */
	    i_arrival_ = prev_pkt_rxTime_ - prev_group_etime_arrival_;
	    i_departure_ = prev_pkt_txTime_ - prev_group_etime_departure_;

	    i_delay_var_ = i_arrival_ - i_departure_;

	    // Save Current group information into prev_group_*.
            prev_group_sseq_ = curr_group_sseq_;
	    prev_group_stime_ = curr_group_stime_;
	    prev_group_eseq_ = prev_pkt_seq_;
	    prev_group_etime_arrival_ = prev_pkt_rxTime_;
	    prev_group_etime_departure_ = prev_pkt_txTime_;

	    group_size_interval_ = prev_group_size_ - curr_group_size_;
	    prev_group_size_ = curr_group_size_;
	    curr_group_size_ = 0;

	    //Update new group information.
	    curr_group_num++;
	    curr_group_sseq_ = sequence;
	    curr_group_stime_ = txTimestamp;

	    m_group_changed_ = true;
	}
    }
    
    curr_group_size_ += packet_size;

    prev_pkt_seq_ = sequence;
    prev_pkt_txTime_ = txTimestamp;
    prev_pkt_rxTime_ = rxTimestamp:;

    return;
}

bool GccRecvController::produceRembFeedback(uint64_t nowUs,
                                      uint16_t sequence,
                                      uint64_t txTimestampMs, uint64_t rxTimestampMs, uint64_t packet_size, uint8_t ecn){
    UpdateGroupInfo(nowUs, sequence, txTimestamp, rxTimestamp, packet_size);
    
    if(m_group_changed_){
        m_group_changed_ = false;
        /**
	 * Produce Delay Based Estimation based on follow variables.
	 * 
	 * i_arrival_ : inter arrival time of two adjacency group.
	 * i_departure_ : inter departure time of two adjacency group.
	 * i_delay_var_ : inter delay variance (inter arrival time - inter departure time) of two adjacency group.
	 * group_size_interval_ : group size interval.
	 *
	 */

	// TODO

    }

}
float GccRecvController::getBandwidth(uint64_t nowUs) const {

    return m_initBw;
}

void GccRecvController::updateMetrics() {
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

void GccRecvController::logStats(uint64_t nowUs) const {

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
