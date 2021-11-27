#pragma once

namespace lipm_walking
{

inline double clamp(double value, double vmin, double vmax){
    if(value > vmax) {
        return vmax;
    }
    if(value < vmin) {
        return vmin;
    }
    return value;
}

}