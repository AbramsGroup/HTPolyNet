# -*- coding: utf-8 -*-
'''
All MD process
'''

import sys
import os
import subprocess
from shutil import copyfile
from shutil import move

class md(object):
    def __init__(self, GMX, MPI, nCPU):
        self.gmx = GMX
        self.mpi = MPI
        self.cpu = nCPU
        
    def checkMDFinish(self, fileName):
        if os.path.isfile('{}.gro'.format(fileName)):
            return True
        else:
            return False
    
    def getLastFrame(self, fileName):
        cmd0 = 'echo 0 > 1.txt'
        subprocess.call(cmd0, shell=True)
        cmd1 = '{} trjconv -s {}.tpr -f {}.trr -o tst.gro -dump -1 < 1.txt'.format(self.gmx, fileName, fileName)
        a1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = a1.communicate()
        a1.wait()
        f1 = open('out1.txt','w')
        f1.write(str(err))
        tmp = 't=0'
        for s in str(err).split():
            if 't=' in s:
                tmp = s
        tmp1 = tmp.strip('t=').split('\\')[0]
        f2 = open('frame.txt','w')
        f3 = open('frame-step.txt', 'a')
        f2.write(tmp1)
        f3.write(cmd1 + '\n')
        f3.write(tmp1 + '\n')
        f2.close()
        f3.close()
        
        cmd2 = '{} trjconv -s {}.tpr -f {}.trr -o {}.gro -dump {} < 1.txt'.format(self.gmx, fileName, fileName, fileName, tmp1)
        a2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        a2.wait()
    
    def extraRun(self, fileName, topName, num):
        self.getLastFrame(fileName)
        name = ''.join([i for i in fileName if not i.isdigit()])
        inName = name + str(num)
        move('{}.gro'.format(fileName), '{}.gro'.format(inName))
    
        outName = name + str(num + 1)
        self.NPTSimulation(inName, topName, outName, 'npt', check=False, re=False)       
        if self.checkMDFinish(outName):
            copyfile('{}.gro'.format(outName), 'npt.gro')
    
    def emSimulation(self, groName, topName, outName, size=True, boxSize=[0, 0, 0], check=True):
        cmd0 = '{} editconf -f {}.gro -box {} {} {} -o {} -c'.format(self.gmx, groName, boxSize[0], boxSize[1], boxSize[2], groName)
        cmd1 = '{} grompp -f em.mdp -c {}.gro -p {}.top -o {} -maxwarn 2'.format(self.gmx, groName, topName, outName)
        cmd2 = '{} -np {} {} mdrun -deffnm {}'.format(self.mpi, self.cpu, self.gmx, outName)
        if size:
            subprocess.call(cmd0, shell=True)
        subprocess.call(cmd1, shell=True)
        subprocess.call(cmd2, shell=True)
        if not self.checkMDFinish(outName):
            cmd3 = '{} -np {} {} mdrun -deffnm {} -rdd 1'.format(self.mpi, self.cpu, self.gmx, outName)
            subprocess.call(cmd3, shell=True)
            if not self.checkMDFinish(outName):
                cmd4 = '{} -np {} {} mdrun -deffnm {} -rdd 0.5'.format(self.mpi, self.cpu, self.gmx, outName)
                subprocess.call(cmd4, shell=True)
                if not self.checkMDFinish(outName):
                    cmd5 = '{} -np {} {} mdrun -deffnm {} -rdd 0.1'.format(self.mpi, self.cpu, self.gmx, outName)
                    subprocess.call(cmd5, shell=True)
                    if not self.checkMDFinish(outName):
                        cmd6 = '{} mdrun -deffnm {}'.format(self.gmx, outName)
                        subprocess.call(cmd6, shell=True)
                        if not self.checkMDFinish(outName):
                            cmd7 = '{} mdrun -deffnm {} -rdd 0.5'.format(self.gmx, outName)
                            subprocess.call(cmd7, shell=True)
                            if not self.checkMDFinish(outName):
                                if check:
                                    sys.exit('Cannot equilibrium well')
                                else:
                                    return False
    
    def NVTSimulation(self, groName, topName, outName, mdpName, check=True):
        cmd1 = '{} grompp -f {}.mdp -c {}.gro -p {}.top -o {} -maxwarn 2'.format(self.gmx, mdpName, groName, topName, outName)
        cmd2 = '{} -np {} {} mdrun -deffnm {}'.format(self.mpi, self.cpu, self.gmx, outName)
        subprocess.call(cmd1, shell=True)
        subprocess.call(cmd2, shell=True)
            
        if not self.checkMDFinish(outName):
            print('NVT didnt finish, please check log file!')
            if check:
                sys.exit()
            else:
                return False
    
    def NPTSimulation(self, groName, topName, outName, mdpName, check=True, re=True):
        cmd1 = '{} grompp -f {}.mdp -c {}.gro -p {}.top -o {} -maxwarn 2'.format(self.gmx, mdpName, groName, topName, outName)
        cmd2 = '{} -np {} {} mdrun -deffnm {}'.format(self.mpi, self.cpu, self.gmx, outName)
        subprocess.call(cmd1, shell=True)
        subprocess.call(cmd2, shell=True)
        
        if not self.checkMDFinish(outName):
            print('NPT didnt finish, please check log file!')
            if not re:
                if check:
                    sys.exit()
                else:
                    return False
            else:
                cmd3 = '{} -np {} {} mdrun -deffnm {} -rdd 1'.format(self.mpi, self.cpu, self.gmx, outName)
                subprocess.call(cmd3, shell=True)
                if not self.checkMDFinish(outName):
                    cmd4 = '{} -np {} {} mdrun -deffnm {} -rdd 0.5'.format(self.mpi, self.cpu, self.gmx, outName)
                    subprocess.call(cmd4, shell=True)
                    if not self.checkMDFinish(outName):
                        cmd5 = '{} -np {} {} mdrun -deffnm {} -rdd 0.1'.format(self.mpi, self.cpu, self.gmx, outName)
                        subprocess.call(cmd5, shell=True)
                        if not self.checkMDFinish(outName):
                            cmd6 = '{} mdrun -deffnm {}'.format(self.gmx, outName)
                            subprocess.call(cmd6, shell=True)
                            if not self.checkMDFinish(outName):
                                cmd7 = '{} mdrun -deffnm {} -rdd 0.5'.format(self.gmx, outName)
                                subprocess.call(cmd7, shell=True)
                                if not self.checkMDFinish(outName):
                                    if check:
                                        sys.exit('Cannot equilibrium well')
                                    else:
                                        return False