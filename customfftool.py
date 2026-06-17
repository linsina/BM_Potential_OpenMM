from openmm import CustomNonbondedForce
from openmm import NonbondedForce
from openmm import Discrete2DFunction
from openmm.app.internal import amber_file_parser
import openmm.unit as units
from math import sqrt




def remove_old_cnbf(system): # Reinitialize is required after this function is called, if it is called after calling "Simulation"(creating context)
     for i in range(system.getNumForces()):
          force = system.getForce(i)
          if isinstance(force, CustomNonbondedForce):
               customnbf = force
               energy_expression = customnbf.getEnergyFunction()
               print(f"\n\nFound old CustomNonbondedForce with energy expression: {energy_expression}\n\n")
               print("Removing old CustomNonbondedForce...\n\n")
               system.removeForce(i)
               break
     if customnbf is None:
          raise RuntimeError("No CustomNonbondedForce found in the system")
     return

def add_126_4_or_126_4_bm_opt(prmtop, system, nonbondedMethod, nonbondedCutoff, use_ig=True, lj_1264=False, lj_1264_bm_opt=True): # Reinitialize is required after this function is called, if it is called after calling "Simulation"(creating context)
     if lj_1264 and lj_1264_bm_opt:
          raise ValueError("Both lj_1264 and bm_opt cannot be True at the same time!")
     elif not lj_1264 and not lj_1264_bm_opt:
          raise ValueError("Both lj_1264 and bm_opt cannot be False at the same time!")
     
     print("WARNING!!!\nWhen using add_126_4_or_126_4_bm_opt make sure the amber atom types of c4 donors(usually metal ions) in your system\ncontain a + or - sign depending on the charge of the donor, otherwise c4 interactions will be calculated wrognly!")
     remove_old_cnbf(system)

     top = amber_file_parser.PrmtopLoader(prmtop)

     ntypes = top.getNumTypes()

     c12_top = [float(x) for x in top._raw_data['LENNARD_JONES_ACOEF']]
     c6_top = [float(x) for x in top._raw_data['LENNARD_JONES_BCOEF']]
     c4_top = [float(x) for x in top._raw_data['LENNARD_JONES_CCOEF']]
     nonb_idx = [int(x) for x in top._raw_data['NONBONDED_PARM_INDEX']]

     c12 = [0 for i in range(ntypes*ntypes)]
     c6 = c12[:]
     c4 = c12[:]

     ene_conv = units.kilocalories_per_mole.conversion_factor_to(units.kilojoules_per_mole)
     length_conv = units.angstroms.conversion_factor_to(units.nanometers)

     c12fac = sqrt(ene_conv) * length_conv**6
     c6fac = ene_conv * length_conv**6
     c4fac = ene_conv * length_conv**4

     if lj_1264_bm_opt:
          a_tot_top = [float(x) for x in top._raw_data['BORN_MAYER_A_TOT']]
          b_top = [float(x) for x in top._raw_data['BORN_MAYER_B']]
          c6_bm_top = [float(x) for x in top._raw_data['BORN_MAYER_C6']]
          c4_bm_top = [float(x) for x in top._raw_data['BORN_MAYER_C4']] 

          a_tot = c12[:]
          b = c12[:]
          c6_bm = c12[:]
          c4_bm = c12[:]

          a_fac = ene_conv
          b_fac = 1/length_conv
          c6_bm_fac = ene_conv * length_conv**6
          c4_bm_fac = ene_conv * length_conv**4

     for i in range(ntypes):
          for j in range(ntypes):
               index = nonb_idx[ntypes*i+j] -1
               if index < 0: continue
               c12[i+ntypes*j] = sqrt(c12_top[index]) * c12fac
               c6[i+ntypes*j] = c6_top[index] * c6fac
               c4[i+ntypes*j] = c4_top[index] * c4fac
               if lj_1264_bm_opt:
                    a_tot[i+ntypes*j] = a_tot_top[index] * a_fac
                    b[i+ntypes*j] = b_top[index] * b_fac
                    c6_bm[i+ntypes*j] = c6_bm_top[index] * c6_bm_fac
                    c4_bm[i+ntypes*j] = c4_bm_top[index] * c4_bm_fac

     ff126 = CustomNonbondedForce('(a/r6)^2 - b/r6; r6=r^6;'
                                   'a=c12(type1, type2);'
                                   'b=c6(type1, type2);')
     ff126.addTabulatedFunction('c12', Discrete2DFunction(ntypes, ntypes, c12))
     ff126.addTabulatedFunction('c6', Discrete2DFunction(ntypes, ntypes, c6))
     ff126.addPerParticleParameter('type')

     ff4 = CustomNonbondedForce('-c/r^4; c=c4(type1, type2)')
     ff4.addTabulatedFunction('c4', Discrete2DFunction(ntypes, ntypes, c4))
     ff4.addPerParticleParameter('type')

     if lj_1264_bm_opt:
          ffbmopt = CustomNonbondedForce('abmtot*exbr - (f6*c6bm)/(r^6) - (f4*c4bm)/(r^4);'
                                             'f6=f4 - exbr*((1/120)*(br^5)*(1 + br/6));'
                                             'f4=1 - exbr*(1 + br*(1 + (1/2)*br*(1 + (1/3)*br*(1 + (1/4)*br))));'
                                             'exbr=exp(-br);'
                                             'br=bbm*r;'
                                             'abmtot=a_tot(type1, type2);'
                                             'bbm=b(type1, type2);'
                                             'c6bm=c6_bm(type1, type2);'
                                             'c4bm=c4_bm(type1, type2);')
          ffbmopt.addTabulatedFunction('a_tot', Discrete2DFunction(ntypes, ntypes, a_tot))
          ffbmopt.addTabulatedFunction('b', Discrete2DFunction(ntypes, ntypes, b))
          ffbmopt.addTabulatedFunction('c6_bm', Discrete2DFunction(ntypes, ntypes, c6_bm))
          ffbmopt.addTabulatedFunction('c4_bm', Discrete2DFunction(ntypes, ntypes, c4_bm))
          ffbmopt.addPerParticleParameter('type')

     for atom in top._getAtomTypeIndexes():
          ff126.addParticle((atom-1,))
          ff4.addParticle((atom-1,))
          if lj_1264_bm_opt:
               ffbmopt.addParticle((atom-1,))
     
     for i in range(system.getNumForces()):
          force = system.getForce(i)
          if isinstance(force, NonbondedForce):
               print("Found NonbondedForce\n\n")
               break

     print("Exceptions are now copied from NonbondedForce to the newly created CustomnoNbondedForces\n\n")
     
     for i in range(force.getNumExceptions()):
          ii, jj, charge, sigma, epsilon = force.getExceptionParameters(i)
          ff126.addExclusion(ii, jj)
          ff4.addExclusion(ii, jj)
          if lj_1264_bm_opt:
               ffbmopt.addExclusion(ii, jj)

     ij_ljtypeidx_lst_prim_ff4 = []
     ij_ljtypeidx_lst_prim_ffbmopt = []

     if use_ig == True:
          for i in range(ntypes):
               for j in range(ntypes):
                    if j >= i:
                         if abs(c4[i+ntypes*j]) > 1e-5:
                              ij_ljtypeidx_lst_prim_ff4.append(i)
                              ij_ljtypeidx_lst_prim_ff4.append(j)
                         if lj_1264_bm_opt:
                              if abs(a_tot[i+ntypes*j]) > 1e-5:
                                   ij_ljtypeidx_lst_prim_ffbmopt.append(i)
                                   ij_ljtypeidx_lst_prim_ffbmopt.append(j)
          
          ij_ljtypeidx_lst_ff4 = set(ij_ljtypeidx_lst_prim_ff4)
          if lj_1264_bm_opt:
               ij_ljtypeidx_lst_ffbmopt = set(ij_ljtypeidx_lst_prim_ffbmopt)

          ig_1_ff4 = []

          for typ in ij_ljtypeidx_lst_ff4:
               for idx, ljtyp in enumerate(top._getAtomTypeIndexes()):
                    if ljtyp-1 == typ:
                         ig_1_ff4.append(idx)

          ig_2_ff4 = []
          ig_1_atypes_ff4 = []
          ig_2_atypes_ff4 = []
          
          to_del_ff4 = []

          for idx, atom in enumerate(ig_1_ff4):
                    atype = top.getAtomTypes()[atom]
                    if "+" in atype or "-" in atype:
                         to_del_ff4.append(idx)
                         ig_2_ff4.append(atom)
                         ig_2_atypes_ff4.append(atype)
                    else:
                         ig_1_atypes_ff4.append(atype)

          for idx in sorted(to_del_ff4, reverse=True):
               del ig_1_ff4[idx]                   

          print("Interaction groups are used for c4 interactions since use_ig = True\n\n")
          print(f"First list in interaction group: {ig_1_atypes_ff4}\n")
          print(f"Second list in interaction group: {ig_2_atypes_ff4}")
          print("The second list should contain only positive and/or negative ion(s), please make sure this is the case!\n\n")

          ff4.addInteractionGroup(ig_1_ff4, ig_2_ff4)

          if lj_1264_bm_opt:
               ig_1_ffbmopt = []

               for typ in ij_ljtypeidx_lst_ffbmopt:
                    for idx, ljtyp in enumerate(top._getAtomTypeIndexes()):
                         if ljtyp-1 == typ:
                              ig_1_ffbmopt.append(idx)

               
               ig_2_ffbmopt = []
               ig_1_atypes_ffbmopt = []
               ig_2_atypes_ffbmopt = []

               to_del_ffbmopt = []

               for idx, atom in enumerate(ig_1_ffbmopt):
                    atype = top.getAtomTypes()[atom]
                    if "+" in atype or "-" in atype:
                         to_del_ffbmopt.append(idx)
                         ig_2_ffbmopt.append(atom)
                         ig_2_atypes_ffbmopt.append(atype)
                    else:
                         ig_1_atypes_ffbmopt.append(atype)

               for idx in sorted(to_del_ffbmopt, reverse=True):
                    del ig_1_ffbmopt[idx]

               print("Interaction groups are used for Born-Mayer Opt force field since use_ig = True\n\n")
               print(f"First list in interaction group: {ig_1_atypes_ffbmopt}\n")
               print(f"Second list in interaction group: {ig_2_atypes_ffbmopt}")
               print("The second list should contain only positive and/or negative ion(s), please make sure this is the case!\n\n")

               ffbmopt.addInteractionGroup(ig_1_ffbmopt, ig_2_ffbmopt)
               
     
     if nonbondedMethod in ('PME', 'LJPME', 'Ewald', 'CutoffPeriodic'):
          ff126.setNonbondedMethod(ff126.CutoffPeriodic)
          ff4.setNonbondedMethod(ff4.CutoffPeriodic)
          ff126.setCutoffDistance(nonbondedCutoff)
          ff4.setCutoffDistance(nonbondedCutoff)
          ff126.setUseLongRangeCorrection(True)
          ff4.setUseLongRangeCorrection(True)
          if lj_1264_bm_opt:
               ffbmopt.setNonbondedMethod(ffbmopt.CutoffPeriodic)
               ffbmopt.setCutoffDistance(nonbondedCutoff)
               ffbmopt.setUseLongRangeCorrection(True)
     elif nonbondedMethod == 'CutoffNonPeriodic':
          ff126.setNonbondedMethod(ff126.CutoffNonPeriodic)
          ff4.setNonbondedMethod(ff4.CutoffNonPeriodic)
          ff126.setCutoffDistance(nonbondedCutoff)
          ff4.setCutoffDistance(nonbondedCutoff)
          if lj_1264_bm_opt:
               ffbmopt.setNonbondedMethod(ffbmopt.CutoffNonPeriodic)
               ffbmopt.setCutoffDistance(nonbondedCutoff)
     elif nonbondedMethod == 'NoCutoff':
          ff126.setNonbondedMethod(ff126.NoCutoff)
          ff4.setNonbondedMethod(ff4.NoCutoff)
          if lj_1264_bm_opt:
               ffbmopt.setNonbondedMethod(ffbmopt.NoCutoff)
     else:
          raise ValueError('Urecognized cutoff option %s' % nonbondedMethod)
     
     system.addForce(ff126)
     system.addForce(ff4)
     if lj_1264_bm_opt:
          system.addForce(ffbmopt)

     return
