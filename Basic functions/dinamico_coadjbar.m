function coadj = dinamico_coadjbar(screw) % optimized on 31.05.2022
coadj = -1*[0 -screw(3) screw(2) 0 -screw(6) screw(5);screw(3) 0 -screw(1) screw(6) 0 -screw(4);-screw(2) screw(1) 0 -screw(5) screw(4) 0;...
        0 -screw(6) screw(5) 0 0 0;screw(6) 0 -screw(4) 0 0 0;-screw(5) screw(4) 0 0 0 0];