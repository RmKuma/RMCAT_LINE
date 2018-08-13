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
 * Dummy controller (CBR) interface for rmcat ns3 module.
 *
 * @version 0.1.1
 * @author Jiantao Fu
 * @author Sergio Mena
 * @author Xiaoqing Zhu
 */

#ifndef GCC_SENDER_CONTROLLER_H
#define GCC_SENDER_CONTROLLER_H

namespace rmcat {

/**
 * Simplistic implementation of a sender-based congestion controller. The
 * algorithm simply returns a constant, hard-coded bandwidth when queried.
 */
class GccSenderController
{
public:
    /** Class constructor */
     GccSenderController();

    /** Class destructor */
    ~GccSenderController();

    /**
     * Set the current bandwidth estimation. This can be useful in test environments
     * to temporarily disrupt the current bandwidth estimation
     *
     * @param [in] newBw Bandwidth estimation to overwrite the current estimation
     */
    // virtual void setCurrentBw(float newBw);

    /**
     * Reset the internal state of the congestion controller
     */
    void reset();

    /**
     * Simplistic implementation of feedback packet processing. It simply
     * prints calculated metrics at regular intervals
     */
    void UpdateLossBasedBitrate(uint64_t nowMs,
                                float plr);
    
    void ApplyDelayBasedBitrate(float DelayBasedEstimateBitrate);
    /**
     * Simplistic implementation of bandwidth getter. It returns a hard-coded
     * bandwidth value in bits per second
     */
    // virtual float getBandwidth(uint64_t nowUs) const;

    float getBitrate();

private:

    // void updateMetrics();
    void logStats(uint64_t nowUs) const;
    void SetMinMaxBitrate(int min_bitrate, int max_bitrate);

    uint64_t m_lastTimeCalcUs;
    bool m_lastTimeCalcValid;
    
    uint64_t m_QdelayUs; /**< estimated queuing delay in microseconds */
    uint32_t m_ploss;  /**< packet loss count within configured window */
    float m_plr;       /**< packet loss ratio within packet history window */

    uint32_t min_configured_bitrate_bps_;
    uint32_t max_configured_bitrate_bps_;

    float m_sending_rate_;  /* Estimated Sending Bps*/
    
    uint32_t min_bitrate_configured_;
    uint32_t max_bitrate_configured_;
    int64_t last_low_bitrate_log_ms_;

    bool m_Bitrate_valid_;
};

}

#endif /* GCC_SENDER_CONTROLLER_H */
