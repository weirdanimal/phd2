/*
 *  star_shape.h - eccentricity of a star from its pixels
 *  PHD Guiding
 *
 *  Copyright (c) 2026 PHD2 Developers
 *  All rights reserved.
 *
 *  This source code is distributed under the following "BSD" license
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *    Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *    Neither the name of openphdguiding.org nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef STAR_SHAPE_H_INCLUDED
#define STAR_SHAPE_H_INCLUDED

#include <cmath>

// Accumulates the intensity-weighted pixels of a star and reports the
// eccentricity of the equivalent ellipse. Kept free of wxWidgets so it can be
// tested on its own.
class StarShape
{
    unsigned int m_n;
    double m_w, m_wx, m_wy, m_wxx, m_wyy, m_wxy;

public:
    // with fewer pixels than this the shape cannot be judged and is reported as round
    enum
    {
        MIN_PIXELS = 5
    };

    StarShape() : m_n(0), m_w(0.), m_wx(0.), m_wy(0.), m_wxx(0.), m_wyy(0.), m_wxy(0.) { }

    // x, y: pixel position relative to any fixed origin (keep it near the star for precision)
    // w: background-subtracted intensity of the pixel
    void Add(double x, double y, double w)
    {
        ++m_n;
        m_w += w;
        m_wx += w * x;
        m_wy += w * y;
        m_wxx += w * x * x;
        m_wyy += w * y * y;
        m_wxy += w * x * y;
    }

    // Eccentricity of the ellipse with the same second moments as the pixels:
    // sqrt(1 - (minor/major)^2). 0 = round, approaching 1 = a line.
    double Eccentricity() const
    {
        if (m_n < MIN_PIXELS || m_w <= 0.0)
            return 0.0;

        // central second moments
        double mx = m_wx / m_w;
        double my = m_wy / m_w;
        double sxx = m_wxx / m_w - mx * mx;
        double syy = m_wyy / m_w - my * my;
        double sxy = m_wxy / m_w - mx * my;

        // eigenvalues of the covariance matrix are the variances along the major and minor axes
        double mean = 0.5 * (sxx + syy);
        double half_diff = 0.5 * (sxx - syy);
        double r = std::sqrt(half_diff * half_diff + sxy * sxy);
        double major = mean + r;
        double minor = mean - r;

        if (major <= 0.0)
            return 0.0;
        if (minor < 0.0)
            minor = 0.0;

        return std::sqrt(1.0 - minor / major);
    }

    // Eccentricity of the object centred near (cx, cy) in a 16-bit image.
    //
    // Uses the pixels within MEASURE_RADIUS of the centre that are more than 3 sigma above the
    // background, weighted by their brightness above the background. The background is measured
    // in a ring just outside that radius. The radius is deliberately larger than the aperture
    // Star::Find uses for HFD, so an extended object such as a galaxy core is not trimmed to a
    // rounder-looking stub.
    enum
    {
        MEASURE_RADIUS = 10,
        RING_WIDTH = 5
    };

    static double MeasureEccentricity(const unsigned short *pixels, int width, int height, double cx, double cy)
    {
        const int A = MEASURE_RADIUS;
        const int B = MEASURE_RADIUS + RING_WIDTH;
        const int A2 = A * A;
        const int B2 = B * B;

        int px = (int) std::floor(cx + 0.5);
        int py = (int) std::floor(cy + 0.5);
        int start_x = px - B < 0 ? 0 : px - B;
        int end_x = px + B > width - 1 ? width - 1 : px + B;
        int start_y = py - B < 0 ? 0 : py - B;
        int end_y = py + B > height - 1 ? height - 1 : py + B;
        if (end_x < start_x || end_y < start_y)
            return 0.0;

        // background mean and sigma in the ring, iteratively clipping outliers at 2 sigma
        double mean_bg = 0.0, sigma_bg = 0.0;
        for (int iter = 0; iter < 9; iter++)
        {
            double sum = 0.0, sum2 = 0.0;
            unsigned int n = 0;
            for (int y = start_y; y <= end_y; y++)
            {
                for (int x = start_x; x <= end_x; x++)
                {
                    int r2 = (x - px) * (x - px) + (y - py) * (y - py);
                    if (r2 <= A2 || r2 > B2)
                        continue;
                    double val = pixels[y * width + x];
                    if (iter > 0 && (val < mean_bg - 2.0 * sigma_bg || val > mean_bg + 2.0 * sigma_bg))
                        continue;
                    sum += val;
                    sum2 += val * val;
                    ++n;
                }
            }
            if (n < 10)
            {
                if (iter == 0)
                    return 0.0; // not enough background to judge
                break;
            }
            double prev_mean = mean_bg;
            mean_bg = sum / n;
            double var = (sum2 - sum * sum / n) / (n - 1);
            sigma_bg = var > 0.0 ? std::sqrt(var) : 0.0;
            if (iter > 0 && std::fabs(mean_bg - prev_mean) < 0.5)
                break;
        }

        double thresh = mean_bg + 3.0 * sigma_bg;

        StarShape shape;
        for (int y = start_y; y <= end_y; y++)
        {
            for (int x = start_x; x <= end_x; x++)
            {
                int dx = x - px, dy = y - py;
                if (dx * dx + dy * dy > A2)
                    continue;
                double val = pixels[y * width + x];
                if (val < thresh)
                    continue;
                shape.Add(dx, dy, val - mean_bg);
            }
        }
        return shape.Eccentricity();
    }
};

#endif // STAR_SHAPE_H_INCLUDED
