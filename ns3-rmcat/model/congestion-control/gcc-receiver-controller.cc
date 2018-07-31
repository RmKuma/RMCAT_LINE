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
 * Gcc Receiver Side Controller (Delay Based Congestion Control) implementation for gcc ns3 module.
 *
 * @version 0.1.0
 */
#include "gcc-receiver-controller.h"
#include <sstream>
#include <cassert>
#define BURST_TIME 5


namespace rmcat {

GccRecvController::GccRecvController() :
    estimated_SendingBps_{0},      // Initialized Estimated Sending Bps TODO 
    
    m_lastTimeCalcUs{0},
    m_lastTimeCalcValid{false},
    m_QdelayUs{0},
    m_ploss{0},
    m_plr{0.f},
    m_RecvR{0.},

    n_loss{0},
    n_total_pkt{0},
    m_lossT{0}, 
    
    valid_pkt_sequence_{0},

    curr_group_num_{0},
    curr_group_sseq_{0},
    curr_group_stime_{0},
    curr_group_eseq_{0},
    curr_group_etime_arrival_{0},
    curr_group_etime_departure_{0},

    prev_group_sseq_{0},
    prev_group_stime_{0},
    prev_group_eseq_{0},
    prev_group_etime_arrival_{0},
    prev_group_etime_departure_{0},

    curr_group_size_{0},
    prev_group_size_{0},
    group_size_interval_{0},

    prev_pkt_seq_{0},
    prev_pkt_txTime_{0},
    prev_pkt_rxTime_{0},

    i_arrival_{0},
    i_departure_{0},
    i_delay_var_{0},

    m_group_changed_{false},
    m_first_packet_{false}{},

    num_of_deltas_{0},
    slope_{8.0/512.0},
    offset_{0},
    prev_offset{0},
    E_{},
    process_offset_{},
    avg_noise_{0.0},
    var_noise_{50},
    ts_delta_hist_{}
    

    {
        E_[0][0] = 100;           
        E_[1][1] = 1e-1;          
        E_[0][1] = E_[1][0] = 0;  
        process_noise_[0] = 1e-13;
        process_noise_[1] = 1e-3; 

    }

GccRecvController::~GccRecvController() {
    ts_delta_hist_.clear();
}

/*void GccRecvController::setCurrentBw(float newBw) {
    m_initBw = newBw;
}
*/

void GccRecvController::reset() {
    m_lastTimeCalcUs = 0;
    m_lastTimeCalcValid = false;

    m_QdelayUs = 0;
    m_ploss = 0;
    m_plr = 0.f;
    m_RecvR = 0.;

    curr_group_num_ = 0;
        
}

float GccRecvController::getBitrate() {
    return estimated_SendingBps_;
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
	t
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

void GccRecvController::UpdateDelayBasedBitrate(uint64_t nowUs,
                                      uint16_t sequence,
                                      uint64_t txTimestampMs, uint64_t rxTimestampMs, uint64_t packet_size, 
				      uint64_t rxRecv_rate, uint8_t ecn){
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
		
		UpdateEstimator(i_arrival_, i_departure_, group_size_interval_, nowUs * 1000);

    }

}

/*
float GccRecvController::getBandwidth(uint64_t nowUs) const {

    return m_initBw;
}
*/

/*
void GccRecvController::updateNetMetrics() {

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

void UpdateEstimator(int t_delta, double ts_delta, int size_delta, int nowMs){
    const double min_frame_period = UpdateMinFramePeriod(ts_delta);
	const double t_ts_delta = t_delta - ts_delta;
	double fs_delta = size_delta;

	++num_of_deltas_;
	if(num_of_deltas_ > 1000){ //1000 is Max value of num_of_deltas
		num_of_deltas_ = 1000;
	}

	//Update Kalman Filter.
    E_[0][0] += process_noise_[0];
    E_[1][1] += process_noise_[1];

	if ((Hypothesis == 'O' &&  offset_ < prev_offset_) ||
    (Hypothesis == 'U' &&  offset_ > prev_offset_)) {
      E_[1][1] += 10 * process_noise_[1];
    }

    const double h[2] = {fs_delta, 1.0};
    const double Eh[2] = {E_[0][0] * h[0] + E_[0][1] * h[1],
                          E_[1][0] * h[0] + E_[1][1] * h[1]};


    const double residual = t_ts_delta - slope_ * h[0] - offset_;
    const bool in_stable_state =  (Hypothesis == 'N');
    const double max_residual = 3.0 * sqrt(var_noise_);
    // We try to filter out very late frames. For instance periodic key
    // frames doesn't fit the Gaussian model well.
    if (fabs(residual) < max_residual) {
       UpdateNoiseEstimate(residual, min_frame_period, in_stable_state);
    } else {
       UpdateNoiseEstimate(residual < 0 ? -max_residual : max_residual,
                      min_frame_period, in_stable_state);
    }

    const double denom = var_noise_ + h[0] * Eh[0] + h[1] * Eh[1];

    const double K[2] = {Eh[0] / denom, Eh[1] / denom};

    const double IKh[2][2] = {{1.0 - K[0] * h[0], -K[0] * h[1]},
                          {-K[1] * h[0], 1.0 - K[1] * h[1]}};
    const double e00 = E_[0][0];
    const double e01 = E_[0][1];

	// Update state.
    E_[0][0] = e00 * IKh[0][0] + E_[1][0] * IKh[0][1];
    E_[0][1] = e01 * IKh[0][0] + E_[1][1] * IKh[0][1];
    E_[1][0] = e00 * IKh[1][0] + E_[1][0] * IKh[1][1];
    E_[1][1] = e01 * IKh[1][0] + E_[1][1] * IKh[1][1];

    // The covariance matrix must be positive semi-definite.
    bool positive_semi_definite =
        E_[0][0] + E_[1][1] >= 0 &&
        E_[0][0] * E_[1][1] - E_[0][1] * E_[1][0] >= 0 && E_[0][0] >= 0;

    slope_ = slope_ + K[0] * residual;
    prev_offset_ = offset_;
    offset_ = offset_ + K[1] * residual;
}

double UpdateMinFramePeriod(double ts_delta){
	
    double m_contrme_period = ts_delta;                           
    if (ts_delta_hist_.size() >= 60) {  
      ts_delta_hist_.pop_front();                                 
    }                                                             
    for (const double old_ts_delta : ts_delta_hist_) {            
      min_frame_period = std::min(old_ts_delta, min_frame_period);
    }                                                             
    ts_delta_hist_.push_back(ts_delta);                           
    return min_frame_period;                                      
}

void UpdateNoiseEstimate(double residual, double ts_delta, bool stable_state){
	if (!stable_state) {
      return;
    }

    // Faster filter during startup to faster adapt to the jitter level
    // of the network. |alpha| is tuned for 30 frames per second, but is scaled
    // according to |ts_delta|.
    double alpha = 0.01;
    if (num_of_deltas_ > 10 * 30) {
      alpha = 0.002;
    }
    // Only update the noise estimate if we're not over-using. |beta| is a
    // function of alpha and the time delta since the previous update.
    const double beta = pow(1 - alpha, ts_delta * 30.0 / 1000.0);
    avg_noise_ = beta * avg_noise_ + (1 - beta) * residual;
    var_noise_ = beta * var_noise_ +
                 (1 - beta) * (avg_noise_ - residual) * (avg_noise_ - residual);
    if (var_noise_ < 1) {
      var_noise_ = 1;
    }
}


}
