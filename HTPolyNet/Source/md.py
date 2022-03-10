# -*- coding: utf-8 -*-
'''
All MD process -- from ming 4/28/21
'''

import sys
import os
import subprocess
from shutil import copyfile
from shutil import move


class md(object):
    def __init__(self, GMX, MPI, nCPU, nGPU=0):
        self.gmx = GMX
        self.gmx_mpi = 'gmx_mpi mdrun'
        self.mpi = MPI
        self.cpu = int(nCPU)
        self.gpu = int(nGPU)

    def combineCMD(self, inOptions):
        opt = ''
        for k, v in inOptions.items():
            tmp = ' -{} {}'.format(k, v)
            opt += tmp
        return opt

    def gmxCMD(self, cmdHeader, inOptions, extraOptions={}, file=False, mpi=True):
        if extraOptions:
            for k, v in extraOptions.items():
                inOptions[k] = v

        if cmdHeader == 'mdrun':
            if self.gpu > 0:
                inOptions['nb'] = 'gpu'
            else:
                inOptions['nb'] = 'cpu'
            opt = self.combineCMD(inOptions)
            cmd = ' '.join(['gmx_mpi', cmdHeader, opt])
            if mpi and self.gpu > 0:
                cmd = '{} -np {} {}'.format(self.mpi, self.gpu, cmd)
            elif mpi and self.cpu > 0:
                cmd = '{} -np {} {}'.format(self.mpi, self.cpu, cmd)
            else:
                cmd = cmd # No parallel calculation
        else:
            opt = self.combineCMD(inOptions)
            cmd = ' '.join(['gmx_mpi', cmdHeader, opt])

        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        junk, output = p.communicate()
        if file:
            i = 0
            fName = '{}{}.log'.format(cmdHeader, i)
            while (os.path.isfile(fName)):
                i += 1
                fName = '{}{}.log'.format(cmdHeader, i)

            with open(fName, 'w') as f:
                f.write(output.decode('utf-8'))

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
        self.NPTSimulation(inName, topName, outName, 'npt-init', check=False, re=False)       
        if self.checkMDFinish(outName):
            copyfile('{}.gro'.format(outName), 'npt.gro')
    
    def emSimulation(self, groName, topName, outName, size=True, boxSize=[0, 0, 0], check=True):
        if size:
            editOptions = {}
            editOptions['f'] = '{}.gro'.format(groName)
            editOptions['box'] = '{} {} {}'.format(boxSize[0], boxSize[1], boxSize[2])
            editOptions['o'] = groName
            self.gmxCMD('editconf', editOptions)

        gromppOptions = {}; mdrunOptions = {}

        gromppOptions['f'] = 'em.mdp'
        gromppOptions['c'] = '{}.gro'.format(groName)
        gromppOptions['p'] = '{}.top'.format(topName)
        gromppOptions['o'] = '{}.tpr'.format(outName)
        gromppOptions['maxwarn'] = '2'
        mdrunOptions['deffnm'] = outName

        self.gmxCMD('grompp', gromppOptions)

        prog = [('rdd', '1', True), ('rdd', '0.5', True), ('rdd', '0.1', True),
                ('ntomp', '6', False), ('rdd', '0.5', False), ('rdd', '0.1', False)]

        self.gmxCMD('mdrun', mdrunOptions, file=True, mpi=True)
        iter = 0
        while not self.checkMDFinish(outName) and iter < len(prog):
            opt = {prog[iter][0]: prog[iter][1]}
            self.gmxCMD('mdrun', mdrunOptions, extraOptions=opt, file=True, mpi=prog[iter][2])
            iter += 1

        if not self.checkMDFinish(outName):
            if check:
                sys.exit('Energy minimization failed. Please check!')
            else:
                return False

    def NPTSimulation(self, groName, topName, outName, mdpName, check=True, re=True):
        gromppOptions = {}; mdrunOptions = {}
        gromppOptions['f'] = '{}.mdp'.format(mdpName)
        gromppOptions['c'] = '{}.gro'.format(groName)
        gromppOptions['p'] = '{}.top'.format(topName)
        gromppOptions['o'] = '{}.tpr'.format(outName)
        gromppOptions['maxwarn'] = '2'
        mdrunOptions['deffnm'] = outName
        mdrunOptions['ntomp'] = '1'

        self.gmxCMD('grompp', gromppOptions)

        prog = [('rdd', '1', True), ('rdd', '0.5', True), ('rdd', '0.1', True),
                 ('rdd', '0.5', False), ('rdd', '0.1', False)]

        self.gmxCMD('mdrun', mdrunOptions, mpi=True)
        if not self.checkMDFinish(outName):
            if not re:
                if check:
                    sys.exit('Cannot equilibrium well')
                else:
                    return False
            else:
                iter = 0
                while not self.checkMDFinish(outName) and iter < len(prog):
                    opt = {prog[iter][0]: prog[iter][1]}
                    self.gmxCMD('mdrun', mdrunOptions, extraOptions=opt, file=True, mpi=prog[iter][2])
                    iter += 1

                if not self.checkMDFinish(outName):
                    if check:
                        sys.exit('NPT failed')
                    else:
                        return False
