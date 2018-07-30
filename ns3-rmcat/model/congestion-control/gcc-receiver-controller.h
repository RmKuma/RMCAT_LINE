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
 * Gcc Receiver Side Controller (Delay Based Congestion Control) interface for gcc ns3 module.
 *
 * @version 0.1.0
 */

#ifndef GCC_RECV_CONTROLLER_H
#define GCC_RECV_CONTROLLER_H

namespace rmcat {

class GccRecvController
{
public:
    /** Class constructor */
    GccRecvController();

    /** Class destructor */
    ~GccRecvController();

    /**
     * Set the current bandwidth estimation. This can be useful in test environments
     * to temporarily disrupt the current bandwidth estimation
     *
     * @param [in] newBw Bandwidth estimation to overwrite the current estimation
     */
    // void setCurrentBw(float newBw);

    /**
     * Reset the internal state of the congestion controller
     */
    void reset();


    /**
     * Function for updating delay based controller's bitrate(Ar).
     * This function will process "Overuse Detecting" and "Overuse Estimating" in internet draft.
     */
    void UpdateDelayBasedBitrate(uint64_t nowUs,
                                 uint16_t sequence,
                                 uint64_t txTimestampMs,
				 uint64_t rxTimestampMs,
				 uint64_t packet_size,
                                 uint64_t rxRecv_rate);
    /**
     * Simplistic implementation of bandwidth getter. It returns a hard-coded
     * bandwidth value in bits per second
     */
    // float getBandwidth(uint64_t nowUs) const;
    
    /**
     * Get Funtion of estimated_SendingBps_
     */
    float GetBitrate();

private:
    // void updateNetMetrics();
    void logStats(uint64_t nowUs) const;

    void UpdateGroupInfo(uint64_t nowUs, uint64_t sequence, uint64_t txTimestamp, uint64_t rxTimestamp, uint64_t packet_size);
    
    float estimated_SendingBps_;  /* Sending rate estimated by delay based controller. */

    uint64_t m_lastTimeCalcUs;
    bool m_lastTimeCalcValid;

    uint64_t m_QdelayUs; /**< estimated queuing delay in microseconds */
    uint32_t m_ploss;  /**< packet loss count within configured window */
    float m_plr;       /**< packet loss ratio within packet history window */
    float m_RecvR;     /**< updated receiving rate in bps */

    int curr_group_num_;

    uint16_t curr_group_sseq_;     /* Current group's first packet sequence*/
    uint64_t curr_group_stime_;    /* Current group's first packet txTimestamp*/
    uint16_t curr_group_eseq_;     /* Current group's last packet sequence.*/
    uint64_t curr_group_etime_arrival_;    /* Current group's last packet rxTimestamp.*/
    uint64_t curr_group_etime_departure_;    /* Current group's last packet txTimestamp.*/

    uint16_t prev_group_sseq_;     /* Previous group's first packet sequence*/
    uint64_t prev_group_stime_;    /* Previous group's first packet txTimestamp.*/
    uint16_t prev_group_eseq_;     /* Previous group's last packet sequence.*/
    uint64_t prev_group_etime_departure_;    /* Previous group's last packet rxTimestamp(arrival time of group's last packet).*/
    uint64_t prev_group_etime_arrival_;      /* Previous group's last packet txTimestamp.*/

    int curr_group_size_;
    int prev_group_size_;
    int group_size_interval_;

    uint16_t prev_pkt_seq_;
    uint64_t prev_pkt_txTime_;
    uint64_t prev_pkt_rxTime_;

    int i_arrival_;
    int i_departure_;
    int i_delay_var_;

    bool m_group_changed_;
    bool m_first_packet_;
};

}

#endif /* DUMMY_CONTROLLER_H */
