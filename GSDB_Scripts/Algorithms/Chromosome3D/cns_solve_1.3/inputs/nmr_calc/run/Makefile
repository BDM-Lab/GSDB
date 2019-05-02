###############################################################################
#									      #
#       Makefile for CNSsolve input run directory			      #
#       Paul Adams 2-9-97						      #
#       copyright Yale University                                             #
#									      #
###############################################################################

.SUFFIXES:
.SUFFIXES: .sed

SHELL=/bin/sh

FIRSTFILES=generate_extended

FILES=anneal \
      accept \
      anneal_cv \
      dg_sa \
      ensemble_cv \
      ensemble \
      rmsd_pm \
      pmrefine

LINKS=../../nmr_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

generate_extended.sed:
	@ sed -e \
	  "s/{===>} structure_file=\".*\"/{===>} structure_file=\"trx.mtf\"/; \
	   s/{===>} output_coor=\".*\"/{===>} output_coor=\"trx_extended.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

anneal.sed:
	@ sed -e \
	  "s/{===>} struct.1=\".*\"/{===>} struct.1=\"il8.mtf\"/; \
	   s/{===>} pdb.in.file.1=\".*\"/{===>} pdb.in.file.1=\"il8.pdb\"/; \
	   s/{===>} flg.print.accept=false/{===>} flg.print.accept=true/; \
	   s/{===>} md.type.cool=\".*\"/{===>} md.type.cool=\"cartesian\"/; \
	   s/{===>} pdb.end.count=[0-9]*/{===>} pdb.end.count=5/; \
	   s/{===>} flg.calc.ave.struct=false/{===>} flg.calc.ave.struct=true/; \
	   s/{===>} md.cool.ss=0.015/{===>} md.cool.ss=0.005/; \
	   s/{===>} md.cool.tmpstp=[0-9]*/{===>} md.cool.tmpstp=25/; \
	   s/{===>} md.cart.flag=true/{===>} md.cart.flag=false/; \
	   s/{===>} nmr.noe.file.1=\".*\"/{===>} nmr.noe.file.1=\"il8_noe.tbl\"/; \
	   s/{===>} nmr.noe.ave.mode.1=\".*\"/{===>} nmr.noe.ave.mode.1=\"cent\"/; \
	   s/{===>} nmr.noe.hbnd.file=\".*\"/{===>} nmr.noe.hbnd.file=\"il8_hbonds.tbl\"/; \
	   s/{===>} nmr.noe.ave.mode.hbnd=\".*\"/{===>} nmr.noe.ave.mode.hbnd=\"cent\"/; \
	   s/{===>} nmr.jcoup.file.1=\".*\"/{===>} nmr.jcoup.file.1=\"\"/; \
	   s/{===>} nmr.oneb.file.1=\".*\"/{===>} nmr.oneb.file.1=\"\"/; \
	   s/{===>} nmr.carb.file.1=\".*\"/{===>} nmr.carb.file.1=\"\"/; \
	   s/{===>} nmr.prot.file.1=\".*\"/{===>} nmr.prot.file.1=\"\"/; \
	   s/{===>} nmr.prot.file.2=\".*\"/{===>} nmr.prot.file.2=\"\"/; \
	   s/{===>} nmr.prot.file.3=\".*\"/{===>} nmr.prot.file.3=\"\"/; \
	   s/{===>} nmr.cdih.file=\".*\"/{===>} nmr.cdih.file=\"il8_dihe.tbl\"/; \
	   s/{===>} pdb.out.name=\".*\"/{===>} pdb.out.name=\"anneal_il8\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} par.2=\"\"/{===>} par.2=\"CNS_TOPPAR:axis.param\"/; \
	   s/{===>} struct.1=\".*\"/{===>} struct.1=\"e1_lr.mtf\"/; \
	   s/{===>} pdb.in.file.1=\".*\"/{===>} pdb.in.file.1=\"e1_lr_ref.pdb\"/; \
	   s/{===>} flg.print.accept=false/{===>} flg.print.accept=true/; \
	   s/{===>} pdb.end.count=[0-9]*/{===>} pdb.end.count=5/; \
	   s/{===>} flg.calc.ave.struct=false/{===>} flg.calc.ave.struct=true/; \
	   s/{===>} md.hot.temp=[0-9]*/{===>} md.hot.temp=60000/; \
	   s/{===>} md.cool.temp=[0-9]*/{===>} md.cool.temp=60000/; \
	   s/{===>} md.cart.temp=[0-9]*/{===>} md.cart.temp=3000/; \
	   s/{===>} md.cart.cdih=[0-9]*/{===>} md.cart.cdih=600/; \
	   s/{===>} md.pow.cdih=[0-9]*/{===>} md.pow.cdih=800/; \
	   s/{===>} nmr.noe.file.1=\".*\"/{===>} nmr.noe.file.1=\"e1_lr_noe.tbl\"/; \
	   s/{===>} nmr.noe.file.2=\"\"/{===>} nmr.noe.file.2=\"e1_lr_noe_hb_his.tbl\"/; \
	   s/{===>} nmr.noe.hbnd.file=\".*\"/{===>} nmr.noe.hbnd.file=\"e1_lr_hbond.tbl\"/; \
	   s/{===>} nmr.jcoup.file.1=\".*\"/{===>} nmr.jcoup.file.1=\"e1_lr_phij_nogly.tbl\"/; \
	   s/{===>} nmr.jcoup.file.2=\"\"/{===>} nmr.jcoup.file.2=\"e1_lr_phij_gly.tbl\"/; \
	   s/{===>} nmr.oneb.file.1=\".*\"/{===>} nmr.oneb.file.1=\"\"/; \
	   s/{===>} nmr.carb.file.1=\".*\"/{===>} nmr.carb.file.1=\"e1_lr_cacb_shifts.tbl\"/; \
	   s/{===>} nmr.prot.file.1=\".*\"/{===>} nmr.prot.file.1=\"\"/; \
	   s/{===>} nmr.prot.file.2=\".*\"/{===>} nmr.prot.file.2=\"\"/; \
	   s/{===>} nmr.prot.file.3=\".*\"/{===>} nmr.prot.file.3=\"\"/; \
	   s/{===>} nmr.dani.file.1=\"\"/{===>} nmr.dani.file.1=\"e1_lr_dani.tbl\"/; \
	   s/{===>} nmr.cdih.file=\".*\"/{===>} nmr.cdih.file=\"e1_lr_dihed.tbl\"/; \
	   s/{===>} pdb.out.name=\".*\"/{===>} pdb.out.name=\"anneal_e1_dani\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)
	@ sed -e \
	  "s/{===>} par.1=\".*\"/{===>} par.1=\"CNS_TOPPAR:dna-rna-allatom-hj-opls.param\"/; \
	   s/{===>} struct.1=\".*\"/{===>} struct.1=\"rna.mtf\"/; \
	   s/{===>} pdb.in.file.1=\".*\"/{===>} pdb.in.file.1=\"rna.pdb\"/; \
	   s/{===>} pdb.atom.select=(.*)/{===>} pdb.atom.select=(name O5' or name C5' or name C4' or name C3' or name O3' or name P)/; \
	   s/{===>} flg.print.accept=false/{===>} flg.print.accept=true/; \
	   s/{===>} pdb.end.count=[0-9]*/{===>} pdb.end.count=5/; \
	   s/{===>} flg.calc.ave.struct=false/{===>} flg.calc.ave.struct=true/; \
	   s/{===>} md.hot.temp=[0-9]*/{===>} md.hot.temp=20000/; \
	   s/{===>} md.hot.step=[0-9]*/{===>} md.hot.step=4000/; \
	   s/{===>} md.hot.cdih=[0-9]*/{===>} md.hot.cdih=5/; \
	   s/{===>} md.cool.temp=[0-9]*/{===>} md.cool.temp=20000/; \
	   s/{===>} md.cool.step=[0-9]*/{===>} md.cool.step=4000/; \
	   s/{===>} md.cart.temp=[0-9]*/{===>} md.cart.temp=1000/; \
	   s/{===>} md.cart.step=[0-9]*/{===>} md.cart.step=6000/; \
	   s/{===>} nmr.noe.file.1=\".*\"/{===>} nmr.noe.file.1=\"rna_noe.tbl\"/; \
	   s/{===>} nmr.noe.file.2=\"\"/{===>} nmr.noe.file.2=\"rna_unoe.tbl\"/; \
	   s/{===>} nmr.noe.hbnd.file=\".*\"/{===>} nmr.noe.hbnd.file=\"rna_hbnd.tbl\"/; \
	   s/{===>} nmr.jcoup.file.1=\".*\"/{===>} nmr.jcoup.file.1=\"\"/; \
	   s/{===>} nmr.oneb.file.1=\".*\"/{===>} nmr.oneb.file.1=\"\"/; \
	   s/{===>} nmr.carb.file.1=\".*\"/{===>} nmr.carb.file.1=\"\"/; \
	   s/{===>} nmr.prot.file.1=\".*\"/{===>} nmr.prot.file.1=\"\"/; \
	   s/{===>} nmr.prot.file.2=\".*\"/{===>} nmr.prot.file.2=\"\"/; \
	   s/{===>} nmr.prot.file.3=\".*\"/{===>} nmr.prot.file.3=\"\"/; \
	   s/{===>} nmr.cdih.file=\".*\"/{===>} nmr.cdih.file=\"rna_dihed.tbl\"/; \
	   s/{===>} nmr.plan.file=\".*\"/{===>} nmr.plan.file=\"rna_plan.tbl\"/; \
	   s/{===>} pdb.out.name=\".*\"/{===>} pdb.out.name=\"anneal_rna\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_3.inp)

accept.sed:

anneal_cv.sed:

dg_sa.sed:
	@ sed -e \
	  "s/{===>} par.1=\".*\"/{===>} par.1=\"CNS_TOPPAR:dna-rna-allatom-hj-opls.param\"/; \
	   s/{===>} struct.1=\".*\"/{===>} struct.1=\"rna.mtf\"/; \
	   s/{===>} pdb.in.file.1=\".*\"/{===>} pdb.in.file.1=\"rna.pdb\"/; \
	   s/{===>} flg.dgsa.flag=true/{===>} flg.dgsa.flag=false/; \
	   s/{===>} pdb.end.count=[0-9]*/{===>} pdb.end.count=5/; \
	   s/{===>} flg.calc.ave.struct=false/{===>} flg.calc.ave.struct=true/; \
	   s/{===>} md.dg.type=\".*\"/{===>} md.dg.type=\"complete\"/; \
	   s/{===>} md.dg.metr.flag=false/{===>} md.dg.metr.flag=true/; \
	   s/{===>} md.dg.ref=\".*\"/{===>} md.dg.ref=\"coordinate\"/; \
	   s/{===>} nmr.noe.file.1=\".*\"/{===>} nmr.noe.file.1=\"rna_noe.tbl\"/; \
	   s/{===>} nmr.noe.file.2=\"\"/{===>} nmr.noe.file.2=\"rna_unoe.tbl\"/; \
	   s/{===>} nmr.noe.hbnd.file=\".*\"/{===>} nmr.noe.hbnd.file=\"rna_hbnd.tbl\"/; \
	   s/{===>} nmr.cdih.file=\".*\"/{===>} nmr.cdih.file=\"rna_dihed.tbl\"/; \
	   s/{===>} nmr.plan.file=\".*\"/{===>} nmr.plan.file=\"rna_plan.tbl\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

ensemble_cv.sed:
	@ sed -e \
	  "s/{===>} pdb.end.count=[0-9]*/{===>} pdb.end.count=30/; \
	   s/{===>} nmr.ens.loop.flag=false/{===>} nmr.ens.loop.flag=true/; \
	   s/{===>} nmr.noe.loop.file.1=\".*\"/{===>} nmr.noe.loop.file.1=\"ambtv_noes.tbl\"/; \
	   s/{===>} nmr.noe.file.1=\".*\"/{===>} nmr.noe.file.1=\"ambtv_other.tbl\"/; \
	   s/{===>} nmr.jcoup.file.1=\".*\"/{===>} nmr.jcoup.file.1=\"ambtv_jcoup_ml\"/; \
	   s/{===>} nmr.jcoup.file.2=\".*\"/{===>} nmr.jcoup.file.2=\"ambtv_jcoupgly_ml\"/; \
	   s/{===>} nmr.cdih.file=\".*\"/{===>} nmr.cdih.file=\"ambtv_dihed_ml\"/; \
	   s/{===>} pdb.out.name=\".*\"/{===>} pdb.out.name=\"ambtv_ml_cv\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

ensemble.sed:
	@ sed -e \
	  "s/{===>} nmr.ens.loop.flag=false/{===>} nmr.ens.loop.flag=true/; \
	   s/{===>} nmr.noe.loop.file.1=\".*\"/{===>} nmr.noe.loop.file.1=\"ambtv_noes.tbl\"/; \
	   s/{===>} nmr.noe.file.1=\".*\"/{===>} nmr.noe.file.1=\"ambtv_other.tbl\"/; \
	   s/{===>} nmr.jcoup.file.1=\".*\"/{===>} nmr.jcoup.file.1=\"ambtv_jcoup_ml\"/; \
	   s/{===>} nmr.jcoup.file.2=\".*\"/{===>} nmr.jcoup.file.2=\"ambtv_jcoupgly_ml\"/; \
	   s/{===>} nmr.cdih.file=\".*\"/{===>} nmr.cdih.file=\"ambtv_dihed_ml\"/; \
	   s/{===>} pdb.out.name=\".*\"/{===>} pdb.out.name=\"ambtv_ml\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

rmsd_pm.sed:
	@ sed -e \
	  "s/{===>} nmr.ens.loop.flag=false/{===>} nmr.ens.loop.flag=true/; \
	   s/{===>} pdb.ens.name=\".*\"/{===>} pdb.ens.name=\"ambtv_ml\"/; \
	   s/{===>} pdb.out.name=\".*\"/{===>} pdb.out.name=\"ambtv_ml_sngl\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

pmrefine.sed:
	@ sed -e \
	  "s/{===>} nmr.ens.loop.flag=false/{===>} nmr.ens.loop.flag=true/; \
	   s/{===>} nmr.noe.loop.file.1=\".*\"/{===>} nmr.noe.loop.file.1=\"ambtv_noes.tbl\"/; \
	   s/{===>} nmr.noe.file.1=\".*\"/{===>} nmr.noe.file.1=\"ambtv_other.tbl\"/; \
	   s/{===>} nmr.jcoup.file.1=\".*\"/{===>} nmr.jcoup.file.1=\"ambtv_jcoup_ml\"/; \
	   s/{===>} nmr.jcoup.file.2=\".*\"/{===>} nmr.jcoup.file.2=\"ambtv_jcoupgly_ml\"/; \
	   s/{===>} nmr.cdih.file=\".*\"/{===>} nmr.cdih.file=\"ambtv_dihed_ml\"/; \
	   s/{===>} pdb.ens.name=\".*\"/{===>} pdb.ens.name=\"ambtv_ml\"/; \
	   s/{===>} pdb.in.name=\".*\"/{===>} pdb.in.name=\"ambtv_ml_sngl\"/; \
	   s/{===>} pdb.out.name=\".*\"/{===>} pdb.out.name=\"ambtv_ml_pm\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
