//
//  main.cpp
//  disp_calc
//
//  Created by Michael Sachs on 10/21/13.
//  Need to make sure c++ standard library is libstdc++
//

#include <iostream>
#include "quakelib/QuakeLib.h"

int main(int argc, const char * argv[])
{

    // insert code here...
    quakelib::Element<4> ele;
    quakelib::Vec<3> pt0(0.0, 0.0, 0.0);
    quakelib::Vec<3> pt1(0.0, 0.0, 1000.0);
    quakelib::Vec<3> pt2(1000.0, 0.0, 1000.0);
    quakelib::Vec<3> pt3(1000.0, 0.0, 0.0);
    
    ele.set_rake(50.0);
    ele.set_vert(0, pt0);
    ele.set_vert(1, pt1);
    ele.set_vert(2, pt2);
    ele.set_vert(3, pt3);
    std::cout << ele.rake() << " " << ele.rake_vector() << std::endl;
    return 0;
}

