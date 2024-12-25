/*
 * =====================================================================================
 *
 *       Filename:  koppen.c
 *
 *    Description:  Compute Koppen classification.
 *
 *        Version:  1.0
 *        Created:  07/09/17 10:33:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno Galdon
 *   Organization:  University of Castilla-La Mancha
 *
 *
 * 05/2021 correction
 * The order of calculations is important...
 * Computing B climate types first (precipitation threshold).
 * Based on Rohli et al. (2015)
 * Amended in 2021 by Andres Navarro
 * =====================================================================================
 */
#include <stdlib.h>
#include <assert.h>

#include "log.h"
#include "utils.h"
#include "struct.h"
#include "koppen.h"

/* *
 * Compute Koppen classes from a sample grid, as in Lohman.
 *  grid: Sample grid, dimensions (nlat, nlon)
 *  nlat: Latitude dimension of the grid
 *  nlon: Longitude dimension of the grid
 * */
int koppen(stats_t *grid, size_t nlat, size_t nlon)
{
    assert(grid);
    assert(nlat > 0);
    assert(nlon > 0);

    float r_d;

    // Iterate grid
    for(int i = 0; i < nlat * nlon; i++)
    {
        // Set default missing class
        grid[i].class = CLASS_MASK;

        // Missing value?
        if(grid[i].t_avg <= MISSING_VAL) continue;

        // Compute r_d threshold
        if(grid[i].r_w > (grid[i].r_annual*0.7f))      r_d = (2.f * grid[i].t_avg);
        else if(grid[i].r_s > (grid[i].r_annual*0.7f)) r_d = (2.f * grid[i].t_avg) + 28.f;
        else                                           r_d = (2.f * grid[i].t_avg) + 14.f;

        // B - Dry Climates
        if(grid[i].r_avg <= r_d)
        {
            if(grid[i].r_avg >= r_d / 2.f) grid[i].class = CLASS_BS;
            else grid[i].class = CLASS_BW;
        }
        // A - Tropical Rainy Climates
        else if(grid[i].t_min >= 18.f)
        {
            if(grid[i].r_min >= 6.f) grid[i].class = CLASS_AF;
            else grid[i].class = CLASS_AW;
        }
        // C - Humid Mesothermal Climates
        else if(grid[i].t_min >= -3.f && grid[i].t_min < 18.f)
        {
            if(grid[i].r_wmax >= 3.f * grid[i].r_smin) grid[i].class = CLASS_CS;
            if(grid[i].r_smax >= 10.f * grid[i].r_wmin) grid[i].class = CLASS_CW;
            else grid[i].class = CLASS_CF;
        }
        // D - Humid Microthermal Climates
        else if(grid[i].t_min < -3.f && grid[i].t_max >= 10.f)
        {
            if(grid[i].r_smax >= 10.f * grid[i].r_wmin) grid[i].class = CLASS_DW;
            else grid[i].class = CLASS_DF;
        }
        // E - Polar Climates
        else if(grid[i].t_max < 10.f)
        {
            if(grid[i].t_max >= 0.f && grid[i].t_max < 10.f) grid[i].class = CLASS_ET;
            else grid[i].class = CLASS_EF;
        }
    }

    return(0);
}

/* *
 * Compute Extended Koppen classes from a sample grid, as in Trewartha.
 *  grid: Sample grid (dimensions nlat * nlon)
 *  nlat: Latitude dimension of the grid
 *  nlon: Longitude dimension of the grid
 * */
int koppen_extended(stats_t *grid, size_t nlat, size_t nlon)
{
    assert(grid);
    assert(nlat > 0);
    assert(nlon > 0);

    float r_d;

    // Iterate grid
    for(int i = 0; i < nlat * nlon; i++)
    {
        // Set default missing class
        grid[i].class = CLASS_MASK;

        // Missing value?
        if(grid[i].t_avg <= MISSING_VAL) continue;

        // Compute r_d threshold
        if(grid[i].r_w > (grid[i].r_annual*0.7f))      r_d = (2.f * grid[i].t_avg);
        else if(grid[i].r_s > (grid[i].r_annual*0.7f)) r_d = (2.f * grid[i].t_avg) + 28.f;
        else                                           r_d = (2.f * grid[i].t_avg) + 14.f;

        // B - Dry Climates
        if(grid[i].r_avg <= r_d)
        {
            if(grid[i].r_avg >= r_d / 2.f)
            {
                // BS
                if(grid[i].t_avg > 18.f) grid[i].class = CLASS_EX_BSh;
                else if(grid[i].r_smax > 3.f * grid[i].r_wmax) grid[i].class = CLASS_EX_BSs;
                else if(grid[i].r_wmax > 3.f * grid[i].r_smax) grid[i].class = CLASS_EX_BSw;
                // TODO: BSn - frequent fog
                else grid[i].class = CLASS_EX_BSk;

            }
            else
            {
                // BW
                if(grid[i].t_avg > 18.f) grid[i].class = CLASS_EX_BWh;
                else if(grid[i].r_smax > 3.f * grid[i].r_wmax) grid[i].class = CLASS_EX_BWs;
                else if(grid[i].r_wmax > 3.f * grid[i].r_smax) grid[i].class = CLASS_EX_BWw;
                // TODO: BWn - frequent fog
                else grid[i].class = CLASS_EX_BWk;
            }
        }
        // A - Tropical rainy climates
        else if(grid[i].t_min >= 18.f)
        {
            if(grid[i].r_min >= 6.f) grid[i].class = CLASS_EX_Af; // Af
            else
            {
                float a = 3.94 - mm_to_inches(grid[i].r_annual)/25.f;
                if(grid[i].r_min > a) grid[i].class = CLASS_EX_Am; // Am
                else
                {
                    if(grid[i].r_smin < 6.f) grid[i].class = CLASS_EX_As; // As
                    else
                    {
                        if(grid[i].t_max - grid[i].t_min < 5.f) grid[i].class = CLASS_EX_Ai; // Ai
                        // TODO: Ag - hottest month comes before the solstice
                        grid[i].class = CLASS_EX_Aw; // Aw
                    }
                }
            }
        }
        // C - Humid Mesothermal Climates
        else if(grid[i].t_min >= -3.f && grid[i].t_min < 18.f && grid[i].t_max > 10.f)
        {
            if(grid[i].r_wmax >= 3.f * grid[i].r_smin)
            {
                // Cs
                if(grid[i].t_max > 22.f) grid[i].class = CLASS_EX_Csa;
                else grid[i].class = CLASS_EX_Csb;
                // TODO: Csc - less than 4 months over 10 degrees
            }
            else if(grid[i].r_smax >= 10.f * grid[i].r_wmin)
            {
                // Cw
                if(grid[i].t_max > 22.f) grid[i].class = CLASS_EX_Cwa;
                else grid[i].class = CLASS_EX_Cwb;
                // TODO: Csc - less than 4 months over 10 degrees
            }
            else
            {
                // Cf
                if(grid[i].t_max > 22.f) grid[i].class = CLASS_EX_Cfa;
                else grid[i].class = CLASS_EX_Cfb;
                // TODO: Csc - less than 4 months over 10 degrees
            }
        }
        // D - Humid Microthermal Climates
        else if(grid[i].t_min < -3.f && grid[i].t_max >= 10.f)
        {
            if(grid[i].r_smax >= 10.f * grid[i].r_wmin)
            {
                if(grid[i].t_min < -38.f) grid[i].class = CLASS_EX_Dwd;
                else grid[i].class = CLASS_EX_Dw;
            }
            else
            {
                if(grid[i].t_min < -38.f) grid[i].class = CLASS_EX_Dfd;
                else grid[i].class = CLASS_EX_Df;
            }
        }
        // E - Polar Climates
        else if(grid[i].t_max < 10.f)
        {
            if(grid[i].t_max >= 0.f && grid[i].t_max < 10.f) grid[i].class = CLASS_EX_ET;
            else grid[i].class = CLASS_EX_EF;
        }
    }

    return(0);
}
