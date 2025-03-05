
This is for providing more detail information.
dadi simulation (with example based on BGI_1.0 and ZJU1 reference)
basically, there are two steps for replicating dadi simulations:
first, repeatly generate input files of dadi by sampling from SNP callling files(obtained by using GenerateSNPfilefromvcffile.py or GenerateSNPfilefromvcftable.py for vcf files or sql snp databases for populations, respectively)
This's been done through the process: using a sql outer-join like function, merge all selected population vcf files into a big snp table, and sample sites from it then produce dadi input file for each sampling.(detail see article methods)
we proveded an example file m_sp_sx_pk_gy_lcb_jd_sm_ytg_kbe_fanyaonly_newdilutetodensity0.01(BGI_1.0) that is produeced from such process and can be used to conduct the results presenting in article through the below step.
second step,
using the following commands to simulated the mallard and domesticated spliting under 'split_mig_1_IM' model(detail see code) with each parameter assigned the 
lower bound,initial value, upper bound for the model and with the lenghth those snp from and with number of repeated times.
command: python2 life/src/NGS/Analysis/usedadiPy2_7/dadicode.py -n mallard14 26 -n quantizpool 24 -f m_sp_sx_pk_gy_lcb_jd_sm_ytg_kbe_fanyaonly_newdilutetodensity0.01 -m split_mig_1_IM -p nuA 2 1e-05 50 -p s 0.5 0 1 -p TA 0.3 1e-06 10 -p TS 0.1 1e-06 5 -p m12 0.5 1e-06 8 -p m21 0.5 1e-06 8 -T m_d -l 480650.4721472733 -b 100 20

This command can be repeated called by the wrapper bootstrapdadisimulation.py but with random initial values for parameters of the model.
wrapper will automately collected all the simulated estimates and all the arrays for the head map pictures(the dadi code is modified,mainly 'plotting.py', to output the heatmap values. src can be found in this repos and build in your mechine) 

-

As we mentioned in recent editon submitted to a journal, we replicate most analyses on newly widely used reference ZJU1 and included more recently published duck WGS data to confirm our results. Anayses procedures were recorded in our anayses engine and corresponding data and scripts can be found in server host throgh our granted path way to visit. Or contact us for the code and more details to replicate our research. Here present some original exploration records that we didn't organized:

###1.  dadi input file based on ZJU1 were in simulation output/basedonZJU1 

	
    This command produce dadi input file and quantized all vcf files of '-q' parameters into one 178(+1) individuals 'quantizpool' population (last columns of both Allele1/Allele2) in the dadi input file. That is 9 domestic vcf represented 11 domestic populations (newdomesticbreedslcwhitelcwhitelcwhite23.indvd.vcf are 8 domestic breeds regarded as one population)
	python bioinfodevelop/analysisAppEntry/GenerateSNPfilefromvcffile.py -t toplevelDuck_ZJU1ref -c configfiles/ZJUchromosomes_all.bed -R ZJU1/bjduckallchr.fna -q sramplouzr8.indvd.vcf 8 -q newdomesticbreedslcwhitelcwhitelcwhite23.indvd.vcf 21 -q shaoxing33.indvd.vcf 26 -v mallard29.indvd.vcf 50 -v spotbilled21.indvd.vcf 38 -q beijing33.indvd.vcf 33 -q domesticpool/smpool.vcf 18 -q domesticpool/campbellpool.vcf 18 -q domesticpool/jdpool.vcf 18 -q domesticpool/cvbellpool.vcf 18 -q domesticpool/gypool.vcf 18 -n 178 optional -o mld50sp38mpl8newdom21sx26bj33with5pool10minAN -a all -d 0.01
	or atcgtools.py VJgendadi with same parameters. Randomly selected 1 SNP per 100kb from the joint sites of all 2+9=11 vcf files that 
    
    The 18669 dilute SNPs are from 49180767 total SNPs from 1188361409 bp. the dilute length should be 442112 bp. (so below some below -l were given a slightly bigger number)
	use model split_mig_1_IM
    def split_mig_1_IM(params,ns,pts):
        nuA,s,TA,TS,m12,m21=params
        xx=dadi.Numerics.default_grid(pts)
        phi=dadi.PhiManip.phi_1D(xx)
        phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
        phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
        nuM=s*nuA
        nuB=(1-s)*nuA
        phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
        fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
        return fs
    python2 /opt/life/src/NGS/Analysis/usedadiPy2_7/dadicode.py -n mallard29 50 -n quantizpool 168 -f mld50sp38mpl8newdom21sx26bj33with5pool10minAN.dilutetodensity0.01 -m split_mig_1_IM -p nuA 2 1e-05 50 -p s 0.5 0 1 -p TA 0.3 1e-06 10 -p TS 0.1 1e-06 5 -p m12 0.5 1e-06 8 -p m21 0.5 1e-06 8 -T m_d -l 480650.4721472733 -b 100 20
    paramname nuA 3.9587142159702635 effective pop size 1756269.1011152226
    paramname s 0.8733708025529396 effective pop size 387467.76621357753
    paramname TA 3.091375752595244 generation 2742955.120284756
    paramname TS 0.12083172966547798 generation 107213.11095894508
    paramname m12 0.003653999128532162 migration rate 4.1181440492456225e-09
    paramname m21 1.6087777970361392 migration rate 1.813130895322419e-06
	run again
    paramname nuA 2.399793328023928 effective pop size 1969067.6243236514
    paramname s 0.8653256617434447 effective pop size 710013.1186873885
    paramname TA 0.8567599205987555 generation 1405969.5906049213
    paramname TS 0.07073203602481341 generation 116073.46392086109
    paramname m12 0.006605745454071318 migration rate 4.02536298689839e-09
    paramname m21 2.196899504253837 migration rate 1.3387312623298885e-06
	
    use  '-n mallard29 50 -n spotbilled21 38' similate divergence between mallard and spot-billed.
    Nref: 843472.3176510853
    paramname nuA 2.8502176946188693 effective pop size 2404079.724690311
    paramname s 0.576410945130569 effective pop size 486186.6758087336
    paramname TA 0.8889120707075663 generation 1499545.4489354726
    paramname TS 0.057781528237729696 generation 97474.239080199
    paramname m12 0.18380995558680158 migration rate 1.0896027749830508e-07
    paramname m21 0.014562059328125043 migration rate 8.63220939406624e-09

    -n spotbilled21 38 -n quantizpool
    Nref: 824096.7965195675
    paramname nuA 2.82737295917464 effective pop size 2330028.998221871
    paramname s 0.9310433361251284 effective pop size 767269.8307216092
    paramname TA 0.4277567604207248 generation 705025.9519046148
    paramname TS 0.30841179226304455 generation 508322.3400256667
    paramname m12 1.0820667683124379 migration rate 6.565167907959129e-07
    paramname m21 7.130817854172712 migration rate 4.326444347489583e-06
    run again
    Nref: 572409.6566201949
    paramname nuA 3.7109807002509605 effective pop size 2124201.188354823
    paramname s 0.9051633131707848 effective pop size 518124.2212772869
    paramname TA 1.781902303594518 generation 2039956.1714625447
    paramname TS 0.09180414583217324 generation 105099.15918420916
    paramname m12 0.0007922270530301423 migration rate 6.920105591054001e-10
    paramname m21 1.9363459696418366 migration rate 1.6913987624484122e-06
    
    correct the -l: dadicode.py -n mallard29 50 -n spotbilled21 40 -f mld50sp38mpl8newdom21sx26bj33with5pool10minAN.add149of7.dilutetodensity0.01 -m split_mig_1_IM -p nuA 2 1e-05 50 -p s 0.5 0 1 -p TA 0.3 1e-06 10 -p TS 0.1 1e-06 5 -p m12 0.5 1e-06 8 -p m21 0.5 1e-06 8 -T m_sp -l 442113 -b 100 20
    paramname nuA 1.560571355109201 effective pop size 1731181.3936679775
    paramname s 0.04784461960462138 effective pop size 53075.250276426355
    paramname TA 0.7123805156741632 generation 1580523.5562079286
    paramname TS 0.007081657767264637 generation 15711.725225348615
    paramname m12 0.6269533647848965 migration rate 2.8258317287006584e-07
    paramname m21 1.7377053939108007 migration rate 7.832262035840837e-07
    
    Beside using ZJU1 as reference, the other difference compare to our previous work on BGI_1.0 reference is that we used additional two (nine actually, as newdomesticbreedslcwhitelcwhitelcwhite23 comprise 8 breeds) domestic populations. results looks like consist with our previous inference that divergence time between domestic lineages and mallard/spot-billed are longer than it between mallard and spot-billed. slightly difffernce on the exactly value may due to model difference(split_nm,split_m were the model presented in paper), data sets diffferences, random, or other errors.
    
    try to use the orignal samples to see whether the additonal domestic samples will have obvious influence on the divergence time. we mannully quntified those samples' allele frequencies into one domestic population,quantiz7pop.
	awk 'BEGIN {OFS = "\t"}{weighted_sum1 = 0;weighted_sum2 = 0;for(i=0;i<8;i++){allele1 =$(4 + i);allele2 = $(17 + i);total = allele1 + allele2;if (total > 0){f1=allele1 /total;f2 = allele2 / total} else {f1 = f2 = 0}if(i >= 1 && i<= 7) {weighted_sum1 += f1 * 21;weighted_sum2 += f2 * 21}}$15 = int(weighted_sum1 + 0.5)"\t"$15;$28 = int(weighted_sum2 + 0.5)"\t"$28;print $0}' mld50sp38mpl8newdom21sx26bj33with5pool10minAN.dilutetodensity0.01 > ...(after some manully editing on title)
	awk 'BEGIN {OFS = "\t"}{for(i=1;i<=7;i++){allele1 = $(7 + i);allele2 = $(20 + i);total=allele1+allele2;if(total>0){f1 = allele1 / total;f2 = allele2/total}else{f1=f2=0};$(7 + i) = f1*149/7;$(20 + i)= f2 *149/7}$15 = $(8) + $(9) + $(10) + $(11) + $(12) + $(13) + $(14); $28 = $(21) + $(22) + $(23) + $(24) + $(25) + $(26) + $(27);print $0}' mld50sp38mpl8newdom21sx26bj33with5pool10minAN.dilutetodensity0.01 > mld50sp38mpl8newdom21sx26bj33with5pool10minAN.add149of7.dilutetodensity0.01
    
    

    
    use build-in model
    def split_mig (params , ns , pts ):
        nu1 ,nu2 ,T,m = params
        xx = Numerics . default_grid (pts)
        phi = PhiManip . phi_1D (xx)
        phi = PhiManip . phi_1D_to_2D (xx , phi )
        phi = Integration . two_pops (phi , xx , T, nu1 , nu2 , m12 =m, m21=m)
        fs = Spectrum . from_phi (phi , ns , (xx ,xx ))
        return fs
	use original model
    
F3 and PCA are based on the same dataset that joint-calling from GATK of all sampes of 12 chromosomes and after plink r2=0.5 prune.(supplemented with non pruned full SNPs sets)

####2.  F3-statistic

qp3Pop version:701
inbreed set NO
Number of triples 1
snps: 3957376

number of blocks for block jackknife: 258
.			Source1     Source2     Target     f_3        std.err     Z	SNPs
 result:	spotbilled  mallard     domestic   0.090572  0.001885    48.049	3956735
###3.  PCA on ZJU1

![all147indvd r2_0 5](https://github.com/user-attachments/assets/68c10588-c28b-4caa-9a40-7ac5085d7ada)

supplemented non r2 pruned results.(currently of 3 chromosomes with variants called) 
![all147indvd direct](https://github.com/user-attachments/assets/270916b7-3ff9-438c-8d8e-efc232bd670f)

unpruned SNP set showS a more fine-resulotion of PCA

All filled pch represent data from our original data of this article(Shaoxing data has pubed in ...bullshell eggs...)

All pch not filled represent ducks from public data

Green are all mallard ducks

Blue are all spotbilled ducks
'*' astral represent D2L, D2B we used for assembly

pub Doms include: 2 JianAn, 3 Longcheng jade green, 3 Mawang, 3 Putian, 3 Taiwan, 3 Sansui, 3 Youxian.
LCW: Liancheng White, MPL : Maple

two chromosome-level genomes are deposit in the NCBI and some additional information under assembly/.
cite: https://doi.org/10.1101/2020.02.03.933069

Our this article last time was reject due to a reviewer claims that our findings have been reported in somewhere else for one reason, which is not true and they published later than ours preprint version. During the past submission history, I feel this issue is still not well accepted by the community that time and back-to-back research and additional further research are reasonable. The improved version is presented here, Welcome to communicate about technology and academia. The source code for this project was intended to be open, and it was indeed made public in the past, having been cloned by many. However, due to some misconduct that has undermined the integrity of academic communication and harmed our interests, access to the source code is now restricted to users known to us, access to the source code is now restricted to user who we know.
