

#include <math.h>

double solve_dp(double weights[], double profits[],
                unsigned short sets[], unsigned short num_items, unsigned short num_sets,
                unsigned short digits, short total_charge, double total_charge_diff,
                unsigned short solution[]) {

    unsigned long blowup = (unsigned long) pow(10.0, (double) digits);

    unsigned short k;
    unsigned short i;
    unsigned short j;
    unsigned long d;
    unsigned short offset = 0;

    unsigned long w[num_items];

    // transform weights to non-negative integers
    double tc = (double) total_charge;
    double max_sum = 0.0;
    for (k = 0; k < num_sets; k ++)
    {
        double wmin = INFINITY;
        double wmax = -INFINITY;
        for (i = 0; i < sets[k]; i++)
        {
            if (weights[i + offset] < wmin)
            {
                wmin = weights[i + offset];
            }
            if (weights[i + offset] > wmax)
            {
                wmax = weights[i + offset];
            }
        }
        tc -= wmin;
        max_sum += wmax - wmin;
        for (i = 0; i < sets[k]; i++)
        {
            w[i + offset] = lround(blowup * (weights[i + offset] - wmin));
        }
        offset += sets[k];
    }

    // lower and upper capacity limits
    long upper_signed = lround(blowup * (tc + total_charge_diff));
    long lower_signed = lround(blowup * fmax(tc - total_charge_diff, 0.0));

    // check if feasible solutions may exist
    long reachable = lround(blowup * max_sum);
    if (upper_signed < 0 || lower_signed > reachable)
    {
        // sum of min weights over all sets is larger than the upper bound
        // or sum of max weights over all sets is smaller than the lower bound
        return -INFINITY;
    }

    // conversion to unsigned long is now safe
    unsigned long upper = (unsigned long) upper_signed;
    unsigned long lower = (unsigned long) lower_signed;

    // init DP and traceback tables
    double dp[upper + 1];
    unsigned short tb[(upper + 1) * num_sets];
    dp[0] = 0;
    for (d = 1; d <= upper; d++)
    {
        dp[d] = -INFINITY;
    }
    for (d = 0; d < (upper + 1) * num_sets; d++)
    {
        tb[d] = 0;
    }

    // DP
    double max_val;
    unsigned long max_idx = 0;
    offset = 0;
    // iterate over all sets
    for (k = 0; k < num_sets; k++)
    {
        // iterate over all capacities
        d = upper + 1;
        do
        {
            d--;
            // find max. profit with capacity i over all items j in set k
            max_val = -INFINITY;
            for (j = 0; j < sets[k]; j++)
            {
                if (d >= w[j + offset] && dp[d - w[j + offset]] + profits[j + offset] > max_val)
                {
                    max_val = dp[d - w[j + offset]] + profits[j + offset];
                    max_idx = j;
                }
            }
            if (max_val > -INFINITY)
            {
                dp[d] = max_val;
                // copy old traceback indices
                for (i = 0; i < k; i++) {
                    tb[(d * num_sets) + i] = tb[((d - w[max_idx + offset]) * num_sets) + i];
                }
                // add new index to traceback
                tb[(d * num_sets) + k] = max_idx;
            }
            else
            {
                dp[d] = -INFINITY;
            }
        }
        while (d > 0l);
        offset += sets[k];
    }

    // find max profit
    max_val = -INFINITY;
    unsigned long max_pos = 0;
    for (d = lower; d < upper+1; d++)
    {
        if (dp[d] > max_val)
        {
            max_val = dp[d];
            max_pos = d;
        }
    }

    // copy solution
    if (max_val > -INFINITY)
    {
        for (k = 0; k < num_sets; k++)
        {
            solution[k] = tb[(max_pos * num_sets) + k];
        }
    }

    return max_val;
}
