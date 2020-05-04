# SdStaf 1.0.6

## Breaking changes
* The function `varSelection` was update. Also `ThinOcc` was improving.


# SdStaf 1.0.3

## Breaking changes
* The function `reduce.env` was modificated. It's improved with parallal methods using `reduce.env(env, transfer, occ_data, mask, parallel=FALSE )`. This method improve the speed of the processing.

* The function `stim.M` was modificated. The radio of buffer zone can be calculate with three methods. `user` if you define the radio value, `Mx.dist` if the radio is the maximum distance between all points and `mean` if the radio is the mean between all points.


# sdStaf 1.0.2

## Breaking changes

* The function `stability` was modificated. New implementation is made with continue values of suitability species distribution. Use `stability(current, proj, continue=TRUE)`
