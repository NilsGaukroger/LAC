nsec = 27;
n = linspace(1,nsec,nsec);
z = (n-1)/nsec.*(Rnew-rotor.radii(1));
twist=interp1(rotor.radii-rotor.radii(1),rad2deg(result.beta),z+rotor.radii(1));
chord=interp1(rotor.radii-rotor.radii(1),result.c,z+rotor.radii(1));
relthicc=interp1(rotor.radii-rotor.radii(1),result.t./result.c*100,z+rotor.radii(1));
twist(1)=twist(2);
chord(1)=chord(2);
relthicc(1)=relthicc(2);
prebent=[0 7.006e-05 %prebent from the DTU10MW htc file
-2.06477e-05 -0.0122119
-0.0072881 -0.0249251
-0.0189235 -0.0273351
-0.0541282 -0.0282163
-0.126633 -0.021321
-0.225666 -0.0128378
-0.288563 -0.00770659
-0.399194 -0.00488317
 -0.576634 -0.0180296
 -0.707136 -0.0501772
 -0.791081 -0.0941228
 -0.837195 -0.14888
 -0.853948 -0.214514
 -0.849367 -0.290618
 -0.79392 -0.462574
 -0.716284 -0.688437
 -0.634358 -0.960017
 -0.553179 -1.28424
 -0.475422 -1.66402
 -0.40318 -2.10743
 -0.330085 -2.6563
 -0.31014 -2.78882
 -0.286719 -2.92517
 -0.255823 -3.06577
 -0.207891 -3.20952
 -0.089894 -3.33685];
prebent = prebent*Rnew/DTU.R;
printinghtc = [n' prebent(:,1) prebent(:,2) z' twist'];
printingae = [z' chord' relthicc'];
disp('Put this in the .htc at nsec 27:')
sprintf('        sec	%d %f %f %f %f;\n',printinghtc')
disp('Put this in the ae.dat:')
sprintf('%f %f %f 1 ;\n',printingae')


