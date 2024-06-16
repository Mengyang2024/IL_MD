import parmed as pmd
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import time

def smi2xyz(smi, name):
    f = open('%s.xyz'%name, 'w')
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    xyz = Chem.MolToXYZBlock(mol)
    print(xyz, file=f)
    print('xyz file generated')

def gauss_opt(name, charge, spin):
    f_input = open('%s.com'%name, 'w')
    f_input.write('%mem=3200MB\n')
    f_input.write('%'+'nproc=8\n')
    f_input.write('%' + ' chk=%s.chk\n'%name)
    f_input.write('#p HF/6-31g(d) opt freq Pop=MK IOp(6/33=2,6/41=10,6/42=17)\n')
    f_input.write('\n')
    f_input.write('%s\n'%name)
    f_input.write('\n')
    f_input.write('%d %d\n'%(charge, spin))
    f_xyz = open('%s.xyz'%name, 'r')
    lines = f_xyz.readlines()
    for line in lines[2:]:
        f_input.write(line)
    f_input.write('\n')
    print('gaussian input file generated')
    os.system('bsub -n 8 "g16 %s.com"' % name)

def check_gauss_log_file(log_file_path):
    f = open(log_file_path, 'r')
    lines = f.readlines()
    if 'Normal' in lines[-1]:
        for line in lines:
            word = line.split()
            try:
                if word[0] == 'Frequencies':
                    for item in word[2:]:
                        if '-' in item:
                            return 'negativa_frequency'
            except:
                continue
        return 'positive_frequency'
    else:
        return 'error'


def run_gauss_opt(smi, name, charge, spin):
    while True:
        if os.path.exists('%s.log'%name):
            log_file_state = check_gauss_log_file('%s.log' % name)
            if log_file_state == 'error':
                os.makedirs('error_in_gauss_log')
                exit()
            elif log_file_state == 'negativa_frequency':
                os.system('obabel -ig16 %s.log -oxyz -O %s.xyz' % (name, name))
                gauss_opt(name, charge, spin)
                while True:
                    if '%s.com' % name in os.popen("bjobs").read():
                        time.sleep(10)
                    else:
                        break
            else:
                break
        else:
            smi2xyz(smi, name)
            gauss_opt(name, charge, spin)
            while True:
                if '%s.com' % name in os.popen("bjobs").read():
                    time.sleep(10)
                else:
                    break

def read_mol2_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def write_mol2_file(file_path, lines):
    with open(file_path, 'w') as file:
        file.writelines(lines)

def modify_charges(lines):
    in_atom_section = False
    atom_index = 0
    
    for i, line in enumerate(lines):
        if line.startswith('@<TRIPOS>ATOM'):
            in_atom_section = True
            continue
        elif line.startswith('@<TRIPOS>'):
            in_atom_section = False
        
        if in_atom_section and line.strip() != '':
            parts = line.split()
            if len(parts) > 8:
                parts[8] = f'{float(parts[8])*0.8:.6f}'
                lines[i] = ' '.join(parts) + '\n'
                atom_index += 1
    
    return lines

def update_mol2_charges(input_file, output_file):
    lines = read_mol2_file(input_file)
    modified_lines = modify_charges(lines)
    write_mol2_file(output_file, modified_lines)

def gen_IL_ff(cation_smi, anion_smi):
    os.mkdir('cation_ff')
    os.chdir('cation_ff')
    run_gauss_opt(cation_smi,'c_1', 1, 1)
    os.system("antechamber -i c_1.log -fi gout -o c_1.mol2 -fo mol2 -at gaff2 -pf yes -c resp")
    update_mol2_charges('c_1.mol2', 'c_1_edit_charge.mol2')
    os.system("sed -i 's/MOL/CAT/' c_1_edit_charge.mol2")
    os.system("sed -i 's/CATECULE/MOLECULE/' c_1_edit_charge.mol2")
    os.system("parmchk2 -i c_1_edit_charge.mol2 -f mol2 -o c_1.frcmod")
    os.system("tleap -f ../input/leap_cation.in")
    parm = pmd.load_file("c_1.prmtop", "c_1.inpcrd")
    parm.save("cation.top", format='gromacs')
    parm.save("cation.pdb")
    os.chdir('..')
    os.mkdir('anion_ff')
    os.chdir('anion_ff')
    run_gauss_opt(anion_smi,'a_1', -1, 1)
    os.system("antechamber -i a_1.log -fi gout -o a_1.mol2 -fo mol2 -at gaff2 -pf yes -c resp")
    update_mol2_charges('a_1.mol2', 'a_1_edit_charge.mol2')
    os.system("sed -i 's/MOL/ANI/' a_1_edit_charge.mol2")
    os.system("sed -i 's/ANIECULE/MOLECULE/' a_1_edit_charge.mol2")
    os.system("parmchk2 -i a_1_edit_charge.mol2 -f mol2 -o a_1.frcmod")
    os.system("tleap -f ../input/leap_anion.in")
    parm = pmd.load_file("a_1.prmtop", "a_1.inpcrd")
    parm.save("anion.top", format='gromacs')
    parm.save("anion.pdb")
    os.chdir('..')

def gen_system():
    os.mkdir('system')
    os.chdir('system')
    os.system('packmol < ../input/packmol.inp')

    with open('../input/cellulose.top', 'r') as file1:
        cellulose_top = file1.readlines()
    with open('../cation_ff/cation.top', 'r') as file1:
        cation_top = file1.readlines()
    with open('../anion_ff/anion.top', 'r') as file1:
        anion_top = file1.readlines()

    with open('cellulose_IL.top', 'w') as file2:
        defaults_part = False
        for line in cellulose_top:
            if line == '[ defaults ]\n':
                defaults_part = True
            if line == '\n':
                defaults_part = False
            if defaults_part:
                file2.write(line)
        file2.write('\n')
        
        cellulose_atomtypes_part = False
        for line in cellulose_top:
            if line == '[ atomtypes ]\n':
                cellulose_atomtypes_part = True
            if line == '\n':
                cellulose_atomtypes_part = False
            if cellulose_atomtypes_part:
                file2.write(line)

        cation_atomtypes_part = False
        for line in cation_top:
            if line == '\n':
                cation_atomtypes_part = False
            if cation_atomtypes_part:
                file2.write(line)
            if line == '; name    at.num    mass    charge ptype  sigma      epsilon\n':
                cation_atomtypes_part = True

        anion_atomtypes_part = False
        for line in anion_top:
            if line == '\n':
                anion_atomtypes_part = False
            if anion_atomtypes_part:
                file2.write(line)
            if line == '; name    at.num    mass    charge ptype  sigma      epsilon\n':
                anion_atomtypes_part = True
        file2.write('\n')

        cellulose_main_part = False
        for line in cellulose_top:
            if line == '[ moleculetype ]\n':
                cellulose_main_part = True
            if line == '[ system ]\n':
                cellulose_main_part = False
            if cellulose_main_part:
                file2.write(line)
        file2.write('\n')

        cation_main_part = False
        for line in cation_top:
            if line == '[ moleculetype ]\n':
                cation_main_part = True
            if line == '[ system ]\n':
                cation_main_part = False
            if cation_main_part:
                file2.write(line)
        file2.write('\n')

        anion_main_part = False
        for line in anion_top:
            if line == '[ moleculetype ]\n':
                anion_main_part = True
            if line == '[ system ]\n':
                anion_main_part = False
            if anion_main_part:
                file2.write(line)
        file2.write('\n')

        cellulsoe_system_part = False
        for line in cellulose_top:
            if line == '[ system ]\n':
                cellulsoe_system_part = True
            if line == '\n':
                cellulsoe_system_part = False
            if cellulsoe_system_part:
                file2.write(line)
        file2.write('\n')

        rest_part = False
        for line in cellulose_top:
            if line == '[ molecules ]\n':
                rest_part = True
            if line == '\n':
                rest_part = False
            if rest_part:
                file2.write(line)
        file2.write('CAT            1000\n')
        file2.write('ANI            1000\n')
        file2.write('\n')
    os.chdir('..')

if __name__ == '__main__':
    gen_IL_ff('CCN1C=C[N+](C)=C1', 'CC([O-])=O')
    gen_system()
