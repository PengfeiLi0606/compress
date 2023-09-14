# compress

This webpage provides the compress code that implements Young's approach to calculate adiabatic contraction of dark matter halos. The code is maintained by Pengfei Li and Jerry Sellwood. For questions or suggestions, please contact PengfeiLi0606@gmail.com.

If you want to call the code with Python, you can use the following lines:

from subprocess import Popen
Process=Popen('./compress %s %s %s %s %s' % (str(name), str(hmass), str(rs), str(ML), str(MLb)), shell=True )
Process.wait()

Each run will generate a ".fit_quality" file, which tells you the comparison between compressed halos and rotation curves. If you don't want to keep them, you can remove them using the following Python code:

outputname = name+'_'+str(format(ML,'.2f'))+'_'+str(format(MLb,'.2f'))+'_'+str(format(round(hmass, 4), '.4f'))+'_'+str(format(round(rs, 4), '.4f'))
Process=Popen( 'rm output/%s.fitquality' %(str(outputname),), shell=True )

If you use this code, please consider citing the following papers:

(1) The first paper applying the compression code to observational data:

https://ui.adsabs.harvard.edu/abs/2005ApJ...634...70S/abstract

(2) The first paper integrating it into rotation curve fits:

https://ui.adsabs.harvard.edu/abs/2022A%26A...665A.143L/abstract
