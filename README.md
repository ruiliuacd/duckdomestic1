
This is for providing more detail information.
dadi simulation (with example)
basically, there are two steps for replicating dadi simulations:
first, repeatly generate input files of dadi by sampling from SNP callling files(obtained by using GenerateSNPfilefromvcffile.py or GenerateSNPfilefromvcftable.py for vcf files or sql snp databases for populations, respectively)
This's been done through the process: using a sql outer-join like function, merge all selected population vcf files into a big snp table, and sample sites from it then produce dadi input file for each sampling.(detail see article methods)
we proveded an example file m_sp_sx_pk_gy_lcb_jd_sm_ytg_kbe_fanyaonly_newdilutetodensity0.01 that is produeced from such process and can be used to conduct the results presenting in article through the below step.
second step,
using the following commands to simulated the mallard and domesticated spliting under 'split_mig_1_IM' model(detail see code) with each parameter assigned the 
lower bound,initial value, upper bound for the model and with the lenghth those snp from and with number of repeated times.
command: python2 life/src/NGS/Analysis/usedadiPy2_7/dadicode.py -n mallard14 26 -n quantizpool 24 -f m_sp_sx_pk_gy_lcb_jd_sm_ytg_kbe_fanyaonly_newdilutetodensity0.01 -m split_mig_1_IM -p nuA 2 1e-05 50 -p s 0.5 0 1 -p TA 0.3 1e-06 10 -p TS 0.1 1e-06 5 -p m12 0.5 1e-06 8 -p m21 0.5 1e-06 8 -T m_d -l 480650.4721472733 -b 100 20

This command can be repeated called by the wrapper bootstrapdadisimulation.py but with random initial values for parameters of the model.
wrapper will automately collected all the simulated estimates and all the arrays for the head map pictures(the dadi code is modified,mainly 'plotting.py', to output the heatmap values. src can be found in this repos and build in your mechine) 



some of the analyzing process were recorded in http://www.atcgorder.com  and could be easy performed in this web platform. or contact us for the code to fully replicate our research (http://www.atcgorder.cn)
![all147indvd direct](https://github.com/user-attachments/assets/270916b7-3ff9-438c-8d8e-efc232bd670f)

supplemented pub duck sequencing data with our data, performed analysis on ZJU1 beijing duck genome reference.(currently of 3 chromosomes with variants called)

All filled pch represent data from our original data of this article(Shaoxing data has pubed in ...bullshell eggs...)

All pch not filled represent ducks from public data

Green are all mallard ducks

Blue are all spotbilled ducks
'*' astral represent D2L, D2B we used for assembly

pub Doms include: 2 JianAn, 3 Longcheng jade green, 3 Mawang, 3 Putian, 3 Taiwan, 3 Sansui, 3 Youxian.
LCW: Liancheng White, MPL : Maple
