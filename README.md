# compress

This webpage provides the compress code that implements Young's approach to calculate adiabatic contraction of dark matter halos. The code is maintained by Pengfei Li and Jerry Sellwood. For questions or suggestions, please contact PengfeiLi0606@gmail.com.

If you want to call the code with Python, you can use the following lines:

from subprocess import Popen<br>
Process=Popen('./compress %s %s %s %s %s' % (str(name), str(hmass), str(rs), str(ML), str(MLb)), shell=True )<br>
Process.wait()

Each run will generate a ".fit_quality" file, which tells you the comparison between compressed halos and rotation curves. If you don't want to keep them, you can remove them using the following Python code:

outputname = name+'_'+str(format(ML,'.2f'))+'_'+str(format(MLb,'.2f'))+'_'+str(format(round(hmass, 4), '.4f'))+'_'+str(format(round(rs, 4), '.4f'))<br>
Process=Popen( 'rm output/%s.fitquality' %(str(outputname),), shell=True )

If you use this code, please consider citing the following papers:<br>
(1) The paper introducing the new approach to derive dark matter halos from observations:<br>
https://ui.adsabs.harvard.edu/abs/2022A%26A...665A.143L/abstract<br>
(2) The paper quantifying the magnitude of baryonic compression with model galaxies:<br>
https://ui.adsabs.harvard.edu/abs/2022ApJ...927..198L/abstract<br>
(3) The paper introducing the concept of adiabatic contraction of dark matter halos due to baryonic infall:<br>
https://ui.adsabs.harvard.edu/abs/2005ApJ...634...70S/abstract
