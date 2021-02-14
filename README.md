# Build-Model

Providing a method to generate ",lt" file to build model for Moltemplate and Packmol software

You need download the following software: Chem3D(windows) Moltemplate(linux or windows) Packmol(linux or windows) 

The workflow:
(1) Draw a single polymer chain in Chem3D, output ".mol" file 
(2) Use this matlab script to genetrate "single.lt", "single.xyz", "system.inp", and "system.lt"
(3) packmol<sytem.inp
    moltempalte.sh -xyz system.xyz sytem.lt
    cleanup_moltemplate.sh
    
    
    
% Generate single.lt file from Chem3D '*mol' file 
%--------------------------------------------------------------------------
%  input parameter:
%                molfile - the name of the data file output by Chem3D;
%                N_atoms - the number of atoms in a polymer chain 
%                N_name -  define the polymer name 
%                N_chains - the number of polymer chains
%  output file:
%             single.lt - for moltemplate
%             single.xyz -for moltemplate
%             system.lt  -for moltemplate
%             system.inp - for packmol to genetrate system.xyz 
% 
%  create Lammps input file (system.data) in the following way in a 
%  linux terminal system:
%        
%            packmol < system.inp # packmol will generates system.xyz file
%            moltemplate.sh -xyz sytem.xyz system.lt
%            cleanup_motemplate.sh   # 'system.data' will be genetrated
%
%                                                    Lijian Song 2021-02-12
%                                                     BUCT_YingLan Lab  
%--------------------------------------------------------------------------


%------------------------Input parameters----------------------------------

molfile=input('1-please write the name of ".mol" file:','s');
N_atom=input('2-please write the number of atoms in a single polymer chain:');
Name=input('3-please write the name of polymer:','s');
N_chains=input('4-please write the number of chains of polymer:') ;

%-----------Get the total number rows of *.mol file---------------------
fid=fopen(molfile);
row=0;
while ~feof(fid)
    tline=fgetl(fid);
    row=row+1;
end
 fclose(fid);

%--------------Read *.mol file information-------------------------------

fid=fopen(molfile);

atom_data=textscan(fid,'%s %s %d %s %f %f %f %f',N_atom,'HeaderLines',7);
atom_id=atom_data{3};
atom_type=atom_data{4};
atom_x=atom_data{5};
atom_y=atom_data{6};
atom_z=atom_data{7};


fid=fopen(molfile);

N_bond=row-7-N_atom-2-3;
Bond_data=textscan(fid,'%s %s %d %d %d %d',N_bond,'HeaderLines',7+N_atom+2);
bond_id=Bond_data{3};
bond_atomid1=Bond_data{5};
bond_atomid2=Bond_data{6};

%-------------------------------------------------------------------------%
%                 Write information to single.lt file                     %
%-------------------------------------------------------------------------%

%  set up a new file
ltfile=fopen('single.lt','w+');

%  opls-aa and lopls-aa force file
fprintf(ltfile,'%s\n \n%s\n\n','import "oplsaa.lt"','import "loplsaa.lt"');
fprintf(ltfile,'%s\t %s\n\n',Name,'inherits OPLSAA {');

%  write init information
fprintf(ltfile,'%s\n',...
    'atom_style full',...
    'units real',...
    'bond_style harmonic',...
    'angle_style harmonic',...
    'dihedral_style opls',...
    'improper_style harmonic',...
    'pair_style lj/cut/coul/long 12.0 12.0',...
    'pair_modify mix geometric shift yes',...
    'special_bonds lj/coul 0.0 0.0 0.5',...
    'dielectric      1.0',...
    'neighbor       2.0 bin',...
    'neigh_modify every 1 delay 0 check yes',...
    'kspace_style pppm 0.0001' );

% read polymer model
fprintf(ltfile,'\n%s\n','read_data system.data');

%write "Data Atoms" information
fprintf(ltfile,'\n%s\n\n %s%d', 'write("Data Atoms") {');
for i=1:N_atom
    fprintf(ltfile,'%s%d\t %s\t %s\t %s\t %f\t %f\t %f\t\n', '$atom:',i,'$mol:.','@atom:    ','0   ',...
        atom_x(i), atom_y(i), atom_z(i));
end
fprintf(ltfile,'%s\n', '}');

%  write bond information
fprintf(ltfile,'\n%s\n', 'write("Data Bond List") {');
for k=1:N_bond
    fprintf(ltfile,'%s%d %s%d %s%d\n', '$bond:',k,'$atom:',...
        bond_atomid1(k),'$atom:',bond_atomid2(k));
    
end
fprintf(ltfile,'\n%s\n %s\n', '}','}');


%-------------------------------------------------------------------------%
%                 Write information to single.xyz file                    %
%-------------------------------------------------------------------------%
xyzfile=fopen('single.xyz','w');
fprintf(xyzfile,'%d\n\n',N_atom);
for j=1:N_atom
    fprintf(xyzfile,'%d %f %f %f\n',atom_id(j),atom_x(j),atom_y(j),atom_z(j));
end

%-------------------------------------------------------------------------%
%                 Write information to system.inp file                    %
%-------------------------------------------------------------------------%
systemxyzfile=fopen('system.inp','w');
fprintf(systemxyzfile,'%s\n',...
'tolerance 2.0',...
'filetype xyz',...
'output system.xyz',...
'structure single.xyz');
fprintf(systemxyzfile,'\n%s %d\n ','number',N_chains);
fprintf(systemxyzfile,'\n%s\n ',...
  'inside box 0.0 0.0 0.0 100 100 100',...
'end structure');

%-------------------------------------------------------------------------%
%                 Write information to system.lt file                    %
%-------------------------------------------------------------------------%
systemltfile=fopen('system.lt','w+');
fprintf(systemltfile,'%s\n\n','import "single.lt"');
fprintf(systemltfile,'\n%s %s %s%d%s\n\n','box = new',Name,'[',N_chains,']');
fprintf(systemltfile,'\n%s',...
'write_once("Data Boundary") {',...
   '0.0  100  xlo xhi',...
   '0.0  100  ylo yhi',...
   '0.0  100  zlo zhi');

fprintf('-------Finished-------------')
fclose(fid);











