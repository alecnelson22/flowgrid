import os, io
import numpy as np
from datetime import *
import getpass
import string
import copy
import math

from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

import vtk

# import pkg_resources  # part of setuptools
# version = pkg_resources.require("ReGrid")[0].version

f2m = 0.3048  # ft to m


class FlowGrid(object):
    def __init__(self):
        self.skip = 0
        self.Prop = {}
        self.out_dir = 'output'

    def __getitem__(self, key):
        return getattr(self, key)

    def set_out_dir(self, dir):
        self.out_dir = dir

    def _check_out(self, subDir):
        if not os.path.exists(os.path.join(self.out_dir, subDir)):
            os.makedirs(os.path.join(self.out_dir, subDir))

    def cell_verts(self, idx):
        """
        Get xyz coordinates for each vertex of a cell
        Grid geometry must be constructed first

        Parameters
        ----------
        idx : array_like of int
            (Size 3) ijk indices of a cell, starting at 1

        Returns
        -------
        coords : array_like of array_like of float
            Sub-arrays contain xyz coordinates of cell vertices
        """
        coords = []
        xyz = [0, 0, 0]
        cell = self.grid.GetCell(idx[0] - 1, idx[1] - 1, idx[2] - 1)
        p_ids = cell.GetPointIds()
        n_ids = pointIds.GetNumberOfIds()
        for n in range(n_ids):
            p = p_ids.GetId(n)
            self.grid.GetPoint(p, xyz)
            coords.append(copy.deepcopy(xyz))
        return coords

    def centroid(self, coords):
        """
        Compute center coordinates of cell

        Parameters
        ----------
        coords : array_like of array_like of float
             (Size 6) Vertex coordinates of cell

        Returns
        -------
        center : array_like of float
            xyz coordinates of center of cell
        """
        return np.mean(coords, axis=0)

    def exportVTK(self, fname):
        """ Saves the SUTRA grid as a VTK file, either a VTKStructuredGrid (.vts)
            or a VTKUnstructuredGrid (.vtu) depending on mesh type.
            fname = the filename it will be saved at, if no extension is given,
            .vts is appended
        """
        filename, ext = os.path.splitext(fname)
        if self.GridType == "vtkStructuredGrid":
            sWrite = vtk.vtkXMLStructuredGridWriter()
            sWrite.SetInputData(self.Grid)
            sWrite.SetFileName(filename + ".vts")
            sWrite.Write()
        elif self.GridType == "vtkUnstructuredGrid":
            sWrite = vtk.vtkXMLUnstructuredGridWriter()
            sWrite.SetInputData(self.Grid)
            sWrite.SetFileName(filename + ".vtu")
            sWrite.Write()
        else:
            print("Grid type is not recognized")

    def printCOORDS(self, f, p, fstr):
        MAXL = 132
        # if self.skip:
        #    self.skip -= 1
        #    return fstr
        for point in p:
            up = " %2.2f" % (point)
            if len(fstr) + len(up) > MAXL:
                f.write(fstr + "\n")
                fstr = " "
            fstr += up
        return fstr

    def printAC(self, f, p, N, fstr):
        MAXL = 132
        if N == 1:
            up = " %i" % (p)
        else:
            up = " %i*%i" % (N, p)
        if len(fstr) + len(up) > MAXL:
            f.write(fstr + "\n")
            fstr = " "
        fstr += up
        return fstr

    def printPROP(self, f, p, N, fstr):
        MAXL = 132
        if N == 1:
            # up = " %1.4e" %(p) # standard notation
            up = " %1.4e" % (p)  # scientific notation
        else:
            up = " %i*%1.4e" % (N, p)
            # up = " %i*%1.4e" %(N,p) # scientific notation
        if len(fstr) + len(up) > MAXL:
            f.write(fstr + "\n")
            fstr = " "
        fstr += up
        return fstr

    def exportTOUGH2(self, fname):
        """Saves the grid as a fixed format TOUGH(2) grid.
        """
        STR = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        self.ne, self.nn, self.nz = np.array(self.Grid.GetDimensions())  # - 1 #
        filename, ext = os.path.splitext(fname)
        if self.GridType == "vtkStructuredGrid":
            with io.open(filename, 'w', newline='\r\n') as f:
                f.write("ELEME")
                # debug
                f.write(
                    """
                    1        10        20        30        40        50        60        70        80
                    |--------|---------|---------|---------|---------|---------|---------|---------|
                    12345678901234567890123456789012345678901234567890123456789012345678901234567890
                    """)

                ii = 0
                for iy in range(self.nn):
                    for ix in range(self.ne):
                        # f.write(str(iy)+str(ix)+"\n")
                        # first base
                        b2 = ii // (len(STR) * len(STR))
                        b1 = (ii - len(STR) * b2) // len(STR)
                        b0 = ii % len(STR)

                        f.write(STR[b2] + STR[b1] + STR[b0] + "\t" + str(ii) + "\n")
                        ii += 1

    def exportECL(self, fname):
        """ Saves the grid as an ECLIPSE grid. For the purposes of ECLIPSE
        """

        # TODO add consistency of dimensions across the inputs
        self.ne, self.nn, self.nz = np.array(self.Grid.GetDimensions()) - 1  # ECLIPSE
        filename, ext = os.path.splitext(fname)
        if self.GridType == "vtkStructuredGrid":
            with io.open(filename + ".GRDECL", 'w', newline='\r\n') as f:
                f.write('-- Generated [\n')
                f.write('-- Format      : ECLIPSE keywords (grid geometry and properties) (ASCII)\n')
                # f.write('-- Exported by : Petrel 2013.7 (64-bit) Schlumberger\n'
                f.write('-- Exported by : ReGrid v.' + version + "\n")
                f.write('-- User name   : ' + getpass.getuser() + "\n")
                f.write('-- Date        : ' + datetime.now().strftime("%A, %B %d %Y %H:%M:%S") + "\n")
                f.write('-- Project     : ' + "ReGrid project\n")
                f.write('-- Grid        : ' + "Description\n")
                f.write('-- Generated ]\n\n')

                f.write('SPECGRID                               -- Generated : ReGrid\n')
                f.write('  %i %i %i 1 F /\n\n' % (self.ne, self.nn, self.nz))
                f.write('COORDSYS                               -- Generated : ReGrid\n')
                f.write('  1 4 /\n\n')  # what is this line?

                f.write('COORD                                  -- Generated : ReGrid\n')
                nz = self.nz
                fstr = str(" ")

                for iy in range(self.nn):
                    for ix in range(self.ne):
                        p0 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(0)
                        fstr = self.printCOORDS(f, p0, fstr)
                        p1 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(4)
                        fstr = self.printCOORDS(f, p1, fstr)
                    # outside edge on far x
                    p2 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(1)
                    fstr = self.printCOORDS(f, p2, fstr)
                    p3 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(5)
                    fstr = self.printCOORDS(f, p3, fstr)
                # outside edge on far y
                for ix in range(self.ne):
                    p8 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(3)
                    fstr = self.printCOORDS(f, p8, fstr)
                    p9 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(7)
                    fstr = self.printCOORDS(f, p9, fstr)
                # outside edge on far northeast
                p14 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(2)
                fstr = self.printCOORDS(f, p14, fstr)
                p15 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(6)
                fstr = self.printCOORDS(f, p15, fstr)
                f.write(fstr)
                fstr = " "
                f.write(" /")
                f.write("\n")
                f.write("\n")

                f.write('ZCORN                                  -- Generated : ReGrid\n')
                for iz in range(self.nz):
                    for iy in range(self.nn):
                        # front face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(0)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(1)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                        # back face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(3)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(2)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                    # bottom layer
                    for iy in range(self.nn):
                        # front face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(4)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(5)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                        # back face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(7)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(6)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                f.write(fstr)
                fstr = " "
                f.write(" /")
                f.write("\n")
                f.write("\n")
                f.write('ACTNUM                                 -- Generated : ReGrid\n')

                c = -999
                N = 0
                for iac in self.ActiveCells.flatten(order='F'):
                    if iac == c:
                        N += 1
                    else:
                        if c != -999:
                            fstr = self.printAC(f, c, N, fstr)
                        c = iac
                        N = 1
                fstr = self.printAC(f, c, N, fstr)
                f.write(fstr)
                f.write(" /")
                f.write("\n")
                f.write("\n")
        else:
            print("Only structured grids can be converted to ECLIPSE files")

    def exportECLPropertyFiles(self, fname):
        """ Convert any point data to cell data
        """

        # Convert point data to cell data for output
        # verifying if this is necessary or if ECLIPSE can use point attributes
        pointConvert = True
        if pointConvert:
            p2c = vtk.vtkPointDataToCellData()
            p2c.SetInputDataObject(self.Grid)
            p2c.PassPointDataOn()
            p2c.Update()
            self.Grid = p2c.GetOutput()

        filename, ext = os.path.splitext(fname)
        for ia in range(self.Grid.GetCellData().GetNumberOfArrays()):
            prop = self.Grid.GetCellData().GetArray(ia).GetName()
            print("exporting prop", prop)
            if self.GridType == "vtkStructuredGrid":
                with io.open(filename + "prop-" + prop.lower() + ".GRDECL", 'w', newline='\r\n') as f:
                    f.write('-- Generated [\n')
                    f.write('-- Format      : ECLIPSE keywords (grid properties) (ASCII)\n')
                    f.write('-- Exported by : ReGrid v.' + version + "\n")
                    f.write('-- User name   : ' + getpass.getuser() + "\n")
                    f.write('-- Date        : ' + datetime.now().strftime("%A, %B %d %Y %H:%M:%S") + "\n")
                    f.write('-- Project     : ' + "ReGrid project\n")
                    f.write('-- Grid        : ' + "Description\n")
                    f.write('-- Unit system : ' + "ECLIPSE-Field\n")
                    f.write('-- Generated ]\n\n')

                    f.write(prop.upper() + '                                 -- Generated : ReGrid\n')
                    f.write('-- Property name in Petrel : ' + prop + '\n')

                    c = -999.9999
                    N = 0
                    ii = 0
                    fstr = " "
                    for iz in range(self.nz):
                        for iy in range(self.nn):
                            for ix in range(self.ne):
                                # iac = round(self.Grid.GetCellData().GetArray(ia).GetTuple1(ii), 4)
                                iac = '{:0.4e}'.format(self.Grid.GetCellData().GetArray(ia).GetTuple1(ii))
                                print(iac)
                                ii += 1
                                if iac == c:
                                    N += 1
                                else:
                                    if c != -999.9999:
                                        fstr = self.printPROP(f, c, N, fstr)
                                    c = eval(iac)
                                    N = 1
                    fstr = self.printPROP(f, c, N, fstr)
                    f.write(fstr)
                    f.write(" /")
                    f.write("\n")


class GRDECL(FlowGrid):
    """
    GRDECL processes Schlumberger ECLIPSE files
    """

    def __init__(self):
        super(GRDECL, self).__init__()

    def loadNodes(self, fname):
        """
            Reads I, J(max), K
                  iterates through I, then decriments J, increments K
                  I = easting
                  J = northing
                  K = depth or elevation?
        """
        with open(fname, "r") as fp:

            # Read in the header
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == "SPECGRID":
                        self.SPECGRID = np.array(fp.readline().split()[0:3], dtype=int)
                    if item[0] == "COORDSYS":
                        self.COORDSYS = fp.readline().split()
                    if item[0] == "COORD":
                        break

            # Read in the coordinates
            self.coords = []
            for line in fp:
                if line.split()[-1] != "/":
                    item = line.split()
                    for c in item:
                        if '*' in c:
                            cc = c.split('*')
                            for i in range(int(cc[0])):
                                self.coords.append(cc[-1])
                        else:
                            self.coords.append(c)
                else:
                    if len(line.split()) > 1:
                        item = line.split()
                        for i in range(len(item) - 1):
                            cc = item[i]
                            if '*' in cc:
                                ccc = cc.split('*')
                                for j in range(int(ccc[0])):
                                    self.coords.append(ccc[-1])
                            else:
                                self.coords.append(c)
                        break
                    else:
                        break

            # Read in ZCORN
            self.zcorn = []
            i = 0
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == "ZCORN":
                        for line in fp:
                            if line.split():
                                if line.split()[-1] != "/":
                                    self.zcorn += line.split()
                                else:
                                    self.zcorn += line.split()[0:-1]
                                    break
                if len(self.zcorn) > 0:
                    break

            # Read in (in)active cells
            self.active = []
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == "ACTNUM":
                        for line in fp:
                            if line.split():
                                if line.split()[-1] != "/":
                                    c = line.split()
                                    if '*' in c:
                                        cc = c.split('*')
                                        for i in range(float(cc[0])):
                                            self.active += cc[-1]
                                    else:
                                        self.active += c
                                else:
                                    self.active += line.split()[0:-1]
                                    break

        self.coords = np.array(self.coords, dtype=float)
        print(self.coords)

        # In Petrel...
        self.ne = self.SPECGRID[0]  # x  i
        self.nn = self.SPECGRID[1]  # y  j
        self.nz = self.SPECGRID[2]  # z  k

        # build grid
        self.buildGrid(plot=False)
        self.buildActiveCells(plot=False)
        self.buildZGrid(plot=False)
        # self.calculateVolumes(plot=False)
        #
        # Convert to VTK
        self.GridType = "vtkStructuredGrid"
        self.Grid = vtk.vtkStructuredGrid()
        self.Grid.SetDimensions(self.ne+1, self.nn+1, self.nz+1)
        vtk_points = vtk.vtkPoints()
        ve = 1.

        for iz in range(self.nz):
            if iz == 0:
                for iy in range(self.nn+1):
                    for ix in range(self.ne+1):
                        vtk_points.InsertNextPoint( self.X0[ix,iy], \
                                                    self.Y0[ix,iy], \
                                               ve * self.ZZT[iz][ix,iy] )
            for iy in range(self.nn+1):
                for ix in range(self.ne+1):
                    vtk_points.InsertNextPoint( self.X0[ix,iy], \
                                                self.Y0[ix,iy], \
                                           ve * self.ZZB[iz][ix,iy] )
        self.Grid.SetPoints(vtk_points)

        # Add in active cells
        ac = vtk.vtkIntArray()
        ac.SetName( "ActiveCells" )
        for iac in self.ActiveCells.flatten( order='F' ):
            ac.InsertNextTuple1( iac )
        self.Grid.GetCellData().AddArray(ac)

    def buildGrid(self, plot=False):
        """
        Topology of COORD mesh, only describes first layer..

                  8--------10-------12-------14
                 /|       /|       /|       /|
                / |      / |      / |      / |
               0--------2--------4--------6  |
               |  9-----|--11----|--13----|--15
               | /      | /      | /      | /
               |/       |/       |/       |/
               1--------3--------5--------7            7  -->   (2*(NE+1))
                                                      15  -->   (2*(NE+1)*(NN+1))
        """

        print("Constructing grid")
        # print("Grid dims", self.ne, self.nn, self.nz)
        # print("Num points", 2*(self.ne+1)*(self.nn+1)*3, len(self.coords))

        # number of edges
        self.ndx = self.ne + 1
        self.ndy = self.nn + 1
        self.ndz = self.nz + 1

        # extract the triplets
        self.points = {}
        self.points["e"] = self.coords[0::3]
        self.points["n"] = self.coords[1::3]
        self.points["z"] = self.coords[2::3]

        print('points e')
        print(self.points["e"])

        # Here are the coordinates
        self.X0 = np.reshape(self.points["e"][0::2] , (self.ndx,self.ndy), order="F")
        self.Y0 = np.reshape(self.points["n"][0::2] , (self.ndx,self.ndy), order="F")
        self.Z0 = np.reshape(self.points["z"][0::2] , (self.ndx,self.ndy), order="F")

        self.X1 = np.reshape(self.points["e"][1::2] , (self.ndx,self.ndy), order="F")
        self.Y1 = np.reshape(self.points["n"][1::2] , (self.ndx,self.ndy), order="F")
        self.Z1 = np.reshape(self.points["z"][1::2] , (self.ndx,self.ndy), order="F")
        #
        # # visualize
        # if plot:
        #     print("plotting")
        #     fig = plt.figure()
        #     ax = fig.add_subplot(111, projection='3d')
        #     ax.plot_wireframe(f2m*self.X0, f2m*self.Y0, f2m*self.Z0, rstride=1, cstride=1)
        #     ax.plot_wireframe(f2m*self.X1, f2m*self.Y1, f2m*self.Z1, rstride=1, cstride=1)
        #     plt.show()

    def buildZGrid(self, plot=False):
        """
            Petrel provides the ZCORN in a truly arcane ordering--it's awful--and really, the programmers
            deserve a special place in hell for doing this. The ordering is as follows, for a given plane:

             29    36  30   37 31    38 32    39 33    40 34    41 35    42
              _______  _______  ______  _______  _______  _______  _______
             /      / /      / /     / /      / /      / /      / /      /|
            /      / /      / /     / /      / /      / /      / /      / |
           00----01 02----03 04----05 06----07 08----09 10----11 12----13 /
            |  A  | |  B   | |   C  | |   D  | |   E  | |  F   | |   G  |/
           14----15 16----17 18----19 20----21 22----23 24----25 26----27


            This pattern is then repeated for each depth layer, it isn't that clear, but my ASCII art skills
            are already sufficiently challenged.

        """

        print("Constructing Z corners")

        # self.zcorn = np.array(self.zcorn, dtype=float)
        # temp = np.zeros( ((self.ne+1)*(self.nn+1)*self.nz) )
        temp = []
        count = 0
        for item in self.zcorn:

            if "*" in item:
                ct = (int)(item.split("*")[0])
                vl = (float)(item.split("*")[1])
                temp += np.tile(vl, ct).tolist()
                count += ct
            else:
                temp += [(float)(item)]
                count += 1

        # layers = np.resize(temp, (8, self.ne*self.nn*self.nz ))
        layers = np.resize(temp, (self.nz * 2, self.ne * self.nn * 4))
        """
        plt.plot(newtemp[0,:])                    # TOP     0    0
        plt.plot(newtemp[1,:])       # SAME --    # BOTTOM  0    1
        #plt.plot(newtemp[2,:])      # SAME --    # TOP     1    2

        plt.plot(newtemp[3,:])       # SAME --    # BOTTOM  1    3
        #plt.plot(newtemp[4,:])      # SAME --    # TOP     2    4

        plt.plot(newtemp[5,:])       # SAME --    # BOTTOM  2    5
        #plt.plot(newtemp[6,:])      # SAME --    # TOP     3    6
        plt.plot(newtemp[7,:])                    # BOTTOM  3    7
        """
        self.ZZT = {}  # zztop ha ha...two year's later this is still funny -TI
        self.ZZB = {}
        for ilay in range(self.nz):
            self.ZZT[ilay] = np.zeros((self.ndx, self.ndy))
            self.ZZB[ilay] = np.zeros((self.ndx, self.ndy))
            iis = 0
            # plt.plot(layers[ilay*2])
            for iin in range(self.nn):
                nears = {}
                fars = {}
                bnears = {}
                bfars = {}
                for iif in range(2):
                    # top
                    nears[iif] = layers[ilay * 2][iis:iis + 2 * self.ne][0::2].tolist()
                    fars[iif] = layers[ilay * 2][iis:iis + 2 * self.ne][1::2].tolist()
                    layers[ilay * 2][iis:iis + 2 * self.ne][0::2] *= 0.  # check
                    layers[ilay * 2][iis:iis + 2 * self.ne][1::2] *= 0.
                    nears[iif].append(fars[iif][-1])
                    fars[iif] = [nears[iif][0]] + fars[iif]
                    # bottom
                    bnears[iif] = layers[ilay * 2 + 1][iis:iis + 2 * self.ne][0::2].tolist()
                    bfars[iif] = layers[ilay * 2 + 1][iis:iis + 2 * self.ne][1::2].tolist()
                    layers[ilay * 2 + 1][iis:iis + 2 * self.ne][0::2] *= 0.
                    layers[ilay * 2 + 1][iis:iis + 2 * self.ne][1::2] *= 0.
                    bnears[iif].append(bfars[iif][-1])
                    bfars[iif] = [bnears[iif][0]] + bfars[iif]
                    #
                    iis += 2 * self.ne

                self.ZZT[ilay][:, iin] = nears[0]
                self.ZZB[ilay][:, iin] = bnears[0]
                # NaN mask for visualizing, but can be sort of a pain to deal with
                # imask = np.nonzero( 1-self.ActiveCells[:,iin,ilay] )
                # self.ZZT[ilay][:,iin][1::][imask] = np.nan
                # self.ZZB[ilay][:,iin][1::][imask] = np.nan
                # if self.ActiveCells[0,iin,ilay] == 0:
                # self.ZZT[ilay][:,iin][0]  = np.nan
                # self.ZZB[ilay][:,iin][0]  = np.nan
                if iin == self.nn - 1:
                    self.ZZT[ilay][:, iin + 1] = fars[1]
                    self.ZZB[ilay][:, iin + 1] = bfars[1]
                    # NaN mask
                    # self.ZZT[ilay][:,iin+1][1::][imask] = np.nan
                    # self.ZZB[ilay][:,iin+1][1::][imask] = np.nan
                    # if self.ActiveCells[0,iin,ilay] == 0:
                    #    self.ZZT[ilay][:,iin+1][0]  = np.nan
                    #    self.ZZB[ilay][:,iin+1][0]  = np.nan

        print("Layers ||", np.linalg.norm(layers), "||")
        # exit()

        # visualize
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            # ax.plot_wireframe( self.X0, self.Y0, self.Z0, rstride=1, cstride=1)

            ax.plot_wireframe(self.X0, self.Y0, self.ZZT[0], rstride=1, cstride=1, color="blue")
            # ax.plot_wireframe( self.X0, self.Y0, self.ZZT[1], rstride=1, cstride=1, color="blue")
            # ax.plot_wireframe( self.X0, self.Y0, self.ZZT[2], rstride=1, cstride=1, color="blue")
            # ax.plot_wireframe( self.X0, self.Y0, self.ZZT[3], rstride=1, cstride=1, color="blue")

            # ax.plot_wireframe( self.X0, self.Y0, self.ZZB[3], rstride=1, cstride=1, color="green")

            plt.gca().set_xlim(np.min(self.X0), np.max(self.X0))
            plt.gca().set_ylim(np.max(self.Y0), np.min(self.Y0))
            # plt.gca().set_zlim( np.max(self.ZZB[3]),  np.min(self.ZZT[0]) )
            plt.gca().set_zlim(5000, 4000)
            plt.savefig("mesh.png")
            plt.show()

    def buildActiveCells(self, plot=False):

        print("Constructing active cells")
        self.ActiveCells = np.zeros((self.ne * self.nn * self.nz), dtype=int)

        count = 0
        for item in self.active:
            if "*" in item:
                ct = (int)(item.split("*")[0])
                vl = (int)(item.split("*")[1])
                self.ActiveCells[count:count + ct] = vl
                count += ct
            else:
                self.ActiveCells[count] = (int)(item)
                count += 1

        self.ActiveCells = np.reshape(self.ActiveCells, (self.ne, self.nn, self.nz), order="F")

        if plot:
            plt.pcolor(self.X0.T, self.Y0.T, self.ActiveCells[:, :, 0].T, edgecolors='w', linewidths=.1)
            plt.xlabel("easting")
            plt.ylabel("northing")
            plt.gca().set_xlim(np.min(self.X0), np.max(self.X0))
            plt.gca().set_ylim(np.max(self.Y0), np.min(self.Y0))
            plt.gca().xaxis.tick_top()
            plt.gca().xaxis.set_label_position("top")
            plt.show()

    def calculateVolumes(self, plot=False):
        # Iterate over cells, assert that we are dealing with parallelpiped, if so
        #             | u1    u2   u3 |
        #    A = det  | v1    v2   v3 |
        #             | w1    w2   w3 |
        # self.Volumes = 10000*np.random.normal(0,1, (self.ne, self.nn, self.nz) )
        self.Volumes = np.zeros((self.ne, self.nn, self.nz))
        for iiz in range(self.nz):
            for iie in range(self.ne):
                for iin in range(self.nn):

                    if self.ActiveCells[iie, iin, iiz]:

                        u = np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZT[iiz][iie, iin])) - \
                            np.array((self.X0[iie + 1, iin], self.Y0[iie + 1, iin], self.ZZT[iiz][iie, iin]))

                        v = np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZT[iiz][iie, iin])) - \
                            np.array((self.X0[iie, iin + 1], self.Y0[iie, iin + 1], self.ZZT[iiz][iie, iin]))

                        w = np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZT[iiz][iie, iin])) - \
                            np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZB[iiz][iie, iin]))
                        if np.any(u != u) or np.any(v != v) or np.any(w != w):
                            print("NAN!", iie, iin, iiz)
                            exit()
                        V = np.linalg.det(np.array((f2m * u, f2m * v, f2m * w)))
                        self.Volumes[iie, iin, iiz] = np.abs(V)  # in m^3

        vr = ((3. / (4. * np.pi)) * self.Volumes) ** (1. / 3.)  # virtual radius, taking into account porosity

        print("Total grid volume: " + str(np.sum(self.Volumes)) + " m^3")

    def read_prop(self, fname, attr_name):
        """ Reads a single property from a file
        """
        print('Reading ' + attr_name + ' input')
        temp = []
        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] != "--":
                        tag = item[0]
                        break

            for line in fp:
                attribute = line.split()
                if attribute:
                    if attribute[0] != "--":
                        if attribute[-1] != "/":
                            for c in attribute:
                                if '*' in c:
                                    cc = c.split('*')
                                    for i in range(int(cc[0])):
                                        temp.append(cc[-1])
                                else:
                                    temp.append(c)
                        else:
                            attribute.pop()
                            for c in attribute:
                                if '*' in c:
                                    cc = c.split('*')
                                    for i in range(int(cc[0])):
                                        temp.append(cc[-1])
                                else:
                                    temp.append(c)
                            break

                        # #attribute = fp.readline().split()[-1]
                        # attribute = fp.readline().split()
                        # print(attribute)
                        # if attribute[0] != "--":
                        #     self.Prop[tag] = attribute
                        # print("loading", attribute)
                        # for line in fp:
                        #     if line.split():
                        #         if line.split()[0] != "--":
                        #             if line.split()[-1] != "/":
                        #                 temp += line.split()
                        #             else:
                        #                 temp += line.split()[0:-1]
                        #                 break
        print(temp)
        data = np.zeros((self.ne * self.nn * self.nz), dtype=float)
        count = 0
        for item in temp:
            if "*" in item:
                ct = (int)(item.split("*")[0])
                vl = (float)(item.split("*")[1])
                data[count:count + ct] = vl
                count += ct
            else:
                data[count] = (float)(item)
                count += 1

        data = np.reshape(data, (self.ne, self.nn, self.nz), order="F")

        # Add to VTK grid
        ac = vtk.vtkDoubleArray()
        ac.SetName(attr_name)
        for iac in data.flatten(order='F'):
            ac.InsertNextTuple1(iac)
        self.Grid.GetCellData().AddArray(ac)

        return data

    def read_outputs(self, fname, prop_strings, toVTK=True, toNumpy=True):
        """
        Reads per-cell properties from .PRT file for all timesteps

        Pack property keywords to read from .PRT into list of lists
        This saves significant time for large grids as only one pass is required through .PRT file
        Inner lists should contain keywords that enable line denoting property section to be uniquely identified

        :param prop_strings: [[ECL prop keyword, subkey1, subkey2, ...], ...]

        Order props as they appear in .PRT file
        """
        print('Reading output properties\n')
        prop = {}
        prop_idx = 0
        for p in prop_strings:
            prop[p[0]] = {}
        with open(fname, "r") as fp:
            t = 0
            II = []
            build = False
            data = np.zeros(self.ne * self.nn * self.nz)
            for line in fp:
                # Find prop keywords
                if not build:
                    if all(e in line for e in prop_strings[prop_idx]):
                        data = np.zeros(self.ne * self.nn * self.nz)
                        # print('Reading output property: ' + prop_strings[prop_idx][0])
                        # print('t = ' + str(t))
                        build = True
                # Read prop data
                else:
                    item = line.split()
                    if len(item) > 0:
                        if 'I=' in line:
                            II = line.split('I=')[1].split()
                            II = list(map(int, II))
                        elif '(*,' in item[0]:
                            idxs = line.split('(')[1].split(')')[0].replace(',', ' ').split()
                            J = int(idxs[1])
                            K = int(idxs[2])
                            vals = line.split(')')[1].split()
                            for c,I in enumerate(II):
                                if '-' in vals[c]:
                                    vals[c] = 0
                                idx = ((self.ne * self.nn) * (K - 1)) + (self.ne * (J - 1)) + (I - 1)
                                data[idx] = vals[c]
                        elif '--' in item[0]:
                            build = False
                            pname = prop_strings[prop_idx][0]
                            prop[pname] = copy.deepcopy(data)
                            if prop_idx < len(prop_strings) - 1:
                                prop_idx += 1
                            # All properties for current t have been read
                            else:
                                print('Exporting grid for t = ' + str(t))
                                ids = []
                                for pn in prop.keys():
                                    data = prop[pn]
                                    if toNumpy:
                                        if t == 0:
                                            self._check_out(pn)
                                        grid_data = np.reshape(data, (self.ne, self.nn, self.nz), order="F")
                                        np.savez_compressed(os.path.join(self.out_dir, pn, pn + '_' + str(t)), grid_data)
                                    if toVTK:
                                        if t == 0:
                                            self._check_out('vtk')
                                        ac = vtk.vtkDoubleArray()
                                        ac.SetName(pn)
                                        for iac in data:
                                            ac.InsertNextTuple1(iac)
                                        id = self.Grid.GetCellData().AddArray(ac)
                                        ids.append(id)
                                if toVTK:
                                    self.exportVTK(os.path.join(self.out_dir, 'vtk', os.path.basename(os.path.splitext(fname_out)[0], str(t))))
                                    for id in ids:
                                        self.Grid.GetCellData().RemoveArray(id)
                                prop_idx = 0
                                t += 1

    def readWellOutput(self, fname, keys):
        wellOutput = {}
        keyOrder = {}
        readNames = False
        build = False
        skip = 0

        for key in keys:
            wellOutput[key] = {}

        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:

                    # read time series values
                    if build:
                        # I don't know why '1' denotes timestep end in .RSM file
                        if item[0] == '1':
                            build = False
                            continue
                        t = item[0]
                        for idx in keyOrder.keys():
                            if t in wellOutput[keyOrder[idx][0]]:
                                wellOutput[keyOrder[idx][0]][t][keyOrder[idx][1]] = float(item[idx])
                            else:
                                if len(keyOrder[idx]) > 1:
                                    wellOutput[keyOrder[idx][0]][t] = {}
                                    wellOutput[keyOrder[idx][0]][t][keyOrder[idx][1]] = float(item[idx])
                                else:
                                    wellOutput[keyOrder[idx][0]][t] = float(item[idx])
                        continue

                    # get names of wells
                    if readNames:
                        for idx in keyOrder.keys():
                            curr = keyOrder[idx]
                            if keyOrder[idx][0][0] == 'W':
                                keyOrder[idx].append(item[idx - skip])
                        next(fp)
                        next(fp)
                        readNames = False
                        build = True
                        skip = 0
                        continue

                    # find line that contains keys
                    if item[0] == 'TIME':
                        keyOrder = {}
                        # if time found, then keys might be on this same line
                        for j,key in enumerate(keys):
                            for i,var in enumerate(item):
                                if key == var:
                                    if i in keyOrder:
                                        keyOrder[i].append(key)
                                    else:
                                        keyOrder[i] = [key]
                                if j == 0 and var[0] != 'W':
                                    # then it is time related or a field variable
                                    skip += 1
                        if len(keyOrder.keys()) > 0:
                            next(fp)
                            readNames = True
        return wellOutput

    # Exports Numpy array of property (can be wells)
    # TODO: inherit from FlowGrid
    def export_prop(self, title, prop):
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        np.savez_compressed(os.path.join(self.out_dir, title), prop)


class SUTRA(FlowGrid):
    """ SUTRA is a USGS flow modelling code.
    """

    def __init__(self):
        super(SUTRA, self).__init__()
        nx, ny, nz = 0, 0, 0

    def loadNodes(self, fname, nx, ny, nz, ve=-1):
        """ Reads in the points of the grid, ususally in a file called nodewise
            fname = nodes file
            nx = number of cells in the easting(x) direction
            ny = number of cells in the northing (y) direction
            nz = number of cells in depth, positive up
            ve = vertical exaggeration, default is 1 (none)
            This method results in the generation of a VtkStructuredGrid
        """
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.ActiveCells = np.ones((self.nx * self.ny * self.nz), dtype=int)

        X = np.loadtxt(fname, comments="#")
        self.points = np.reshape(np.array((X[:, 2], X[:, 3], X[:, 4])).T, (nx, ny, nz, 3))

        # work directly with VTK structures
        self.GridType = "vtkStructuredGrid"
        self.Grid = vtk.vtkStructuredGrid()
        self.Grid.SetDimensions(nx, ny, nz)
        vtk_points = vtk.vtkPoints()
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    vtk_points.InsertNextPoint(self.points[ix, iy, iz][0], \
                                               self.points[ix, iy, iz][1], \
                                               ve * self.points[ix, iy, iz][2])
        self.Grid.SetPoints(vtk_points)

    def loadNodesConnections(self, nodes, connections):
        """ In contrast to the above method, the points and connections can be loaded instead.
            For non-regular grids this is necessary. This method results in the generation
            of a vtkUnstructuredGrid.
            nodes = node file, often called nodewise
            connections = element connections, often called incident
        """
        X = np.loadtxt(nodes, comments="#")
        #                   x       y       z
        points = np.array((X[:, 2], X[:, 3], X[:, 4])).T

        self.GridType = "vtkUnstructuredGrid"
        self.Grid = vtk.vtkUnstructuredGrid()
        vtk_points = vtk.vtkPoints()

        for point in range(np.shape(points)[0]):
            vtk_points.InsertNextPoint(points[point, 0], points[point, 1], points[point, 2])
        self.Grid.SetPoints(vtk_points)

        # Read in the connections, the format is as follows
        #  nodeid    p0, p1, p2, p3, p4, p5, p6, p7, p8
        C = np.loadtxt(connections, comments="#", skiprows=2, dtype=int)
        for line in range(np.shape(C)[0]):
            idList = vtk.vtkIdList()
            for node in C[line, :][1:]:
                idList.InsertNextId(node - 1)
            self.Grid.InsertNextCell(vtk.VTK_HEXAHEDRON, idList)

    def readPermeability(self, fname, label=("$\kappa_x$", "$\kappa_y$", "$\kappa_z$")):
        """ Reads in SUTRA permeability data
        """
        k = np.loadtxt(fname, comments="#")
        nr, nc = np.shape(k)
        if self.GridType == "vtkStructuredGrid":
            # Sutra and VTK use opposite ordering
            k = np.reshape(k, (self.nx - 1, self.ny - 1, self.nz - 1, np.shape(k)[1]))
            k = np.reshape(k, (nr, nc), order='F')
        kx = vtk.vtkDoubleArray()
        kx.SetName(label[0])
        ky = vtk.vtkDoubleArray()
        ky.SetName(label[1])
        kz = vtk.vtkDoubleArray()
        kz.SetName(label[2])
        for ik, K in enumerate(k):
            kx.InsertNextTuple1(K[2])
            ky.InsertNextTuple1(K[3])
            kz.InsertNextTuple1(K[4])
        self.Grid.GetCellData().AddArray(kx)
        self.Grid.GetCellData().AddArray(ky)
        self.Grid.GetCellData().AddArray(kz)

    def readPorosity(self, fname, label="phi"):  # LaTeX tags work too: $\phi$
        phi = np.loadtxt(fname)
        nr, nc = np.shape(phi)
        if self.GridType == "vtkStructuredGrid":
            # Sutra and VTK use opposite ordering
            phi = np.reshape(phi, (self.nx, self.ny, self.nz, np.shape(phi)[1]))
            phi = np.reshape(phi, (nr, nc), order='F')
        vphi = vtk.vtkDoubleArray()
        vphi.SetName(label)
        for ik, K in enumerate(phi):
            vphi.InsertNextTuple1(K[5])
        self.Grid.GetPointData().AddArray(vphi)

    def readPressure(self, fname, ts=2, label="$P$"):
        nnodes = self.nx * self.ny * self.nz
        P = np.loadtxt(fname, comments="#")[ts * nnodes:(ts + 1) * nnodes, :]
        C = np.loadtxt(fname, comments="#")[ts * nnodes:(ts + 1) * nnodes, :]
        nr, nc = np.shape(P)
        if self.GridType == "vtkStructuredGrid":
            # Sutra and VTK use opposite ordering
            P = np.reshape(P, (self.nx, self.ny, self.nz, np.shape(P)[1]))
            P = np.reshape(P, (nr, nc), order='F')
            C = np.reshape(C, (self.nx, self.ny, self.nz, np.shape(C)[1]))
            C = np.reshape(C, (nr, nc), order='F')
        vP = vtk.vtkDoubleArray()
        vP.SetName(label)

        vC = vtk.vtkDoubleArray()
        vC.SetName("Concentration")

        for ik in range(nnodes):
            vP.InsertNextTuple1(P[ik, 3])
            vC.InsertNextTuple1(C[ik, 4])
            # vP.InsertNextTuple1( P[2*nnodes+ik, 3] )a
        self.Grid.GetPointData().AddArray(vP)
        self.Grid.GetPointData().AddArray(vC)


class CMG(FlowGrid):
    """
    For handling CMG files
    """
    def __init__(self):
        super(CMG, self).__init__()
        self.out_dir = 'output'
        self.out_props = {}

    def CORNER(self, fname, cp_type):
        """
        Builds corner point grid (*GRID *CORNER) geometry from a CMG file

        Parameters
        ----------
        fname : str
            Input file
        cp_type : array_like of str
            *CORNER subkeyword(s), describing how to read corner points
                User should refer to input file to determine appropriate keywords
                Available options are ['CORNERS'] or ['ZCORN', subkey_x, subkey_y]
                    subkey_x options are 'DI' or 'XCORN' (TODO: Implement XCORN)
                    subkey_y options are 'DJ' or 'YCORN' (TODO: Implement YCORN)

            Support has not been added for 'COORD'
        """
        print('Building corner point grid')
        with open(fname, "r") as fp:
            # Read header
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    # Searches for line of format *GRID *CORNER I J K
                    if "GRID" in item[0]:
                        self.gridType = item[1]
                        self.size = np.array(item[2:5], dtype=int)
                        break

            if cp_type[0] == 'CORNERS':
                X, Y, Z = self.read_CORNERS(fp)
                X, Y, Z = self._calc_coords(X, Y, Z)

            elif cp_type[0] == 'ZCORN':
                if cp_type[1] == 'DI':
                    i_widths = self.read_DIDJ(fp, 'I')
                    X = self._write_X(i_widths)
                elif cp_type[1] == 'XCORN':  # TODO: Implement XCORN
                    pass
                if cp_type[2] == 'DJ':
                    j_widths = self.read_DIDJ(fp, 'J')
                    Y = self._write_Y(j_widths)
                elif cp_type[2] == 'YCORN':  # TODO: Implement YCORN
                    pass

                for line in fp:
                    item = line.split()
                    if len(item) > 0:
                        if item[0] == "ZCORN" or item[0] == "*ZCORN":
                            break
                Z = self.read_ZCORN(fp)
                X, Y, Z = self._calc_coords(X, Y, Z)
        self.structured_grid(X, Y, Z)

        # # Read NULL
        # for line in fp:
        #     item = line.split()
        #     # Assumes NULL keyword followed by ALL
        #     if item[0] == "NULL" or item[0] == "*NULL":
        #         break
        # self.buildActiveCells(fp)

        # # Add in active cells
        # ac = vtk.vtkIntArray()
        # ac.SetName("ActiveCells")
        # for iac in self.ActiveCells.flatten(order='F'):
        #     ac.InsertNextTuple1(iac)
        # self.Grid.GetCellData().AddArray(ac)

    # This is lightly tested
    def CART(self, fname):
        """
        Build cartesian grid (*GRID *CART) from CMG file

        Parameters
        ----------
        fname : str
            Input file
        """
        print('Building cartesian grid')
        self.iWidths = []
        self.jWidths = []
        with open(fname, "r") as fp:

            # Read header
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    # Searches for line of format *GRID *CART I J K
                    if item[0] == "GRID" or item[0] == "*GRID":
                        self.gridType = item[1]
                        self.size = np.array(item[2:5], dtype=int)
                        break

            k_spacing = 0
            # Assumes DEPTH is the final keyword describing grid structure (move 'break' if not)
            for line in fp:
                item = [ i.strip('*') for i in line.split() ]
                # Read Z-axis orientation
                if item[0] == "KDIR" and item[1] == "DOWN": # K=1 is top layer
                    kdir = 1
                elif item[0] == "KDIR" and item[1] == "UP": # K=1 is bottom layer
                    kdir = -1
                # Read X-axis spacing
                elif item[0] == "DI":
                    if item[1] == "CON":
                        self.X = np.linspace(0, float(item[2]) * self.size[0], self.size[0]+1)
                # Read Y-axis spacing
                elif item[0] == "DJ":
                    if item[1] == "CON":
                        self.Y = np.linspace(0, float(item[2]) * (self.size[1]), self.size[1]+1)
                # Read Z-axis spacing
                elif item[0] == "DK":
                    if item[1] == "CON":
                        k_spacing = float(item[2])
                # Read DEPTH (assumes of form *DEPTH *TOP I J K depth)
                elif item[0] == "DEPTH":
                    depth = float(item[5])
                    self.Z = np.linspace(depth, depth + (k_spacing * kdir * float(self.size[2])), self.size[2]+1)
                    break

        # Write cell vertex coordinates
        XX, YY, ZZ = ([] for el in range(3))
        for k in range(self.size[2] + 1):
            for j in range(self.size[1] + 1):
                XX.extend(self.X)
                YY.extend([self.Y[j]] * (self.size[0] + 1))
            ZZ.extend([self.Z[k]] * (self.size[0] + 1) * (self.size[1] + 1))
        self.structured_grid(XX, YY, ZZ)

    def read_CORNERS(self, fp):
        """
        Read corner point grid data from keyword *CORNERS

        Parameters
        ----------
        fp : textIOWrapper
            File object to read lines from
            ex) with open(fname, "r") as fp:
                    # move file pointer
                    ...
                    X, Y, Z = read_CORNERS(fp)

        Returns
        -------
        X : arary_like
            x coordinates of corner points
        Y : array_like
            y coordinates of corner points
        Z : array_like
            z coordinates of corner points
        """
        corners = []
        count = 0
        nb = self.size[0] * self.size[1] * self.size[2] * 24

        for line in fp:
            item = line.split()
            if len(item) > 0:
                if "CORNERS" in item[0]:
                    break
        for line in fp:
            item = line.split()
            if len(item) > 0:
                for c in item:
                    if "*" in c:
                        item = c.split("*")
                        for i in range(int(item[0])):
                            corners.append(item[1])
                            count += 1
                    else:
                        corners.append(c)
                        count += 1
                # all attributes have been read
                if count == nb:
                    break
        sp = int(len(corners) / 3)
        X = list(map(float, corners[:sp]))
        Y = list(map(float, corners[sp:2*sp]))
        Z = list(map(float, corners[2*sp:]))
        return X, Y, Z

    def read_DIDJ(self, fp, d):
        """
        Read grid data from keyword *DI or *DJ

        Parameters
        ----------
        fp : textIOWrapper
            File object to read lines from
            ex) with open(fname, "r") as fp:
                    # move file pointer
                    ...
                    read_DIDJ(fp, 'I')
        d : str
            Should be 'I' or 'J', to read either *DI or *DJ
        """
        widths = []
        count = 0
        if d == "I":
            nb = 0
        else:
            nb = 1
        for line in fp:
            item = line.split()
            if len(item) > 0:
                if "D" + d in item[0]:
                    if d + "VAR" in item[1]:
                        widths += item[2:]
                        count += len(item) - 2
                        break
                    elif item[1] == "CON":
                        for i in range(self.size[nb]):
                            widths += [item[2]]
                        return widths
        for line in fp:
            item = line.split()
            if len(item) > 0:
                for zz in item:
                    if "*" in zz:
                        item = zz.split("*")
                        for i in range(0, int(item[0])):
                            widths.append(item[1])
                            count += 1
                    else:
                        widths.append(zz)
                        count += 1
                # If true, all attributes have been read
                if count == self.size[nb]:
                    break
        return widths

    def read_ZCORN(self, fp):
        """
        Read grid data from keyword *ZCORN

        Parameters
        ----------
        fp : textIOWrapper
            File object to read lines from
            ex) with open(fname, "r") as fp:
                    # move file pointer
                    ...
                    Z = read_ZCORN(fp)
        """
        Z = []
        for line in fp:
            item = line.split()
            if len(item) > 0:
                if item[0][0] != "*":
                    for zz in item:
                        if "*" in zz:
                            item = zz.split("*")
                            for i in range(0, int(item[0])):
                                Z.append(float(item[1]))
                        else:
                            Z.append(float(zz))
                else:
                    break
        return Z

    def _write_X(self, i_widths):
        # Converts i_widths to X coords
        X = []
        for layer in range(2):
            for k in range(self.size[2]):
                for j in range(self.size[1]):
                    # NW-T and NE-T corners
                    x1 = 0
                    x2 = 0
                    for i in range(self.size[0]):
                        x2 = x1 + float(i_widths[i])
                        X.extend([x1, x2])
                        x1 = x2
                    # SW-T and SE-T corners
                    x1 = 0
                    x2 = 0
                    for i in range(self.size[0]):
                        x2 = x1 + float(i_widths[i])
                        X.extend([x1, x2])
                        x1 = x2
        return X

    def _write_Y(self, j_widths):
        # Converts j_widths to Y coords
        Y = []
        for layer in range(2):
            for k in range(self.size[2]):
                y = 0
                for j in range(self.size[1]):
                    # NW-T and NE-T corners
                    for i in range(self.size[0]):
                        Y.extend([y, y])
                    y += float(j_widths[j])
                    # SW-T and SE-T corners
                    for i in range(self.size[0]):
                        Y.extend([y, y])
        return Y

    def _calc_coords(self, X, Y, Z):
        """
        Transforms coordinate data following Eclipse ZCORN ordering
        to coordinate data following VTK structured grid ordering
        (see 'Detailed Description' in vtkStructuredGrid class)

        Parameters
        ----------
        Data in input arrays has following structure:
        for nk
            for nj
                NW-T, NE-T, ... for ni
                SW-T, SE-T, ... for ni
            for nj
                NW-B, NE-B, ... for ni
                SW-B, SE-B, ... for ni

        X : array_like
            Cell vertex x coordinates
        Y : array_like
            Cell vertex y coordinates
        Z : array_like
            Cell vertex z coordinates
        """
        def _write_coords(coord):
            XX.append(X[coord])
            YY.append(Y[coord])
            ZZ.append(Z[coord])

        def _build_layer():
            for j in range(self.size[1]):
                for i in range(self.size[0]):
                    # write NW corner
                    if i == 0:
                        nwCoord = 2 * i + 4 * self.size[0] * j + const
                        _write_coords(nwCoord)
                    # write NE corner
                    neCoord = 2 * i + 4 * self.size[0] * j + const + 1
                    _write_coords(neCoord)
                if j == self.size[1] - 1:
                    for i in range(self.size[0]):
                        # write SW corner
                        if i == 0:
                            swCoord = 2 * i + 4 * self.size[0] * j + 2 * self.size[0] + const
                            _write_coords(swCoord)
                        # write SE corner
                        seCoord = 2 * i + 4 * self.size[0] * j + 2 * self.size[0] + const + 1
                        _write_coords(seCoord)

        # At this point, we have all points needed for unstructured grid in X,Y,Z
        # However, they must be re-arranged so we can define Hexahedrons
        # TODO: REFINE CELLS
        # PSUEDO:
        #    find cell to be refined
        #    add new cells (as easy as pie)

        XX, YY, ZZ = ([] for i in range(3))
        const = 0
        for k in range(self.size[2]):
            _build_layer()
            if k == self.size[2] - 1:
                const += self.size[0] * self.size[1] * 4
                _build_layer()
                break
            else:
                const += self.size[0] * self.size[1] * 8
        return XX, YY, ZZ

    def structured_grid(self, X, Y, Z):
        """
        Builds a VTS structured grid
        Assumes X, Y, Z coordinate arrays in correct order
        (see 'Detailed Description' in vtkStructuredGrid class for ordering)

        Parameters
        ----------
        X : array_like
            Cell vertex x coordinates
        Y : array_like
            Cell vertex y coordinates
        Z : array_like
            Cell vertex z coordinates
        """
        self.GridType = "vtkStructuredGrid"
        self.Grid = vtk.vtkStructuredGrid()
        self.Grid.SetDimensions(self.size[0] + 1, self.size[1] + 1, self.size[2] + 1)
        vtk_points = vtk.vtkPoints()
        for point in range(len(X)):
            vtk_points.InsertNextPoint(X[point], Y[point], Z[point])
        self.Grid.SetPoints(vtk_points)

    # ExPeRiMeNtAl - doesn't work yet (but will be pretty darn cool when it does)
    def unstructured_grid(self):
        """
        **************
        IN DEVELOPMENT
        **************

        Transforms a structured grid to an unstructured one
        Unstructured grid will (eventually) allow us to implement
        locally refined grid (LGR) and irregular geometry
        """
        # TODO: change hardcoded 'POR'
        t = vtk.vtkThreshold()
        t.SetInputDataObject(self.Grid)
        t.SetInputArrayToProcess(0, 0, 0, self.Grid.FIELD_ASSOCIATION_CELLS, "POR")
        t.ThresholdByUpper(-1)
        t.Update()

        ug = t.GetOutput()

        test_ids = vtk.vtkIdList()
        test_ids.InsertNextId(1)

        ug.GetCellPoints(1, test_ids)
        ########

        self.GridType = "vtkUnstructuredGrid"
        self.Grid = vtk.vtkUnstructuredGrid()
        vtk_points = vtk.vtkPoints()

        for point in range(np.shape(points)[0]):
            vtk_points.InsertNextPoint(points[point, 0], points[point, 1], points[point, 2])
        self.Grid.SetPoints(vtk_points)

        # Read in the connections, the format is as follows
        #  nodeid    p0, p1, p2, p3, p4, p5, p6, p7, p8
        C = np.loadtxt(connections, comments="#", skiprows=2, dtype=int)
        for line in range(np.shape(C)[0]):
            idList = vtk.vtkIdList()
            for node in C[line, :][1:]:
                idList.InsertNextId(node - 1)
            self.Grid.InsertNextCell(vtk.VTK_HEXAHEDRON, idList)

        ########
        ug = vtk.vtkXMLUnstructuredGridWriter()
        ug.SetDataModeToAscii()
        ug.SetFileName('unstructured_grid.vtu')
        ug.SetInputConnection(t.GetOutputPort())
        ug.Write()

    def buildActiveCells(self, fp):
        """
        fp : textIOWrapper
            File object to read lines from
        """
        self.ActiveCells = []
        count = 0
        for line in fp:
            item = line.split()
            for zz in item:
                if "*" in zz:
                    item = zz.split("*")
                    for i in range(0, int(item[0])):
                        self.ActiveCells.append(int(item[1]))
                        count += 1
                else:
                    self.ActiveCells.append(int(zz))
                    count += 1
            if count == self.size[0] * self.size[1] * self.size[2]:
                break
        self.ActiveCells = np.array(self.ActiveCells)
        self.ActiveCells = np.reshape(self.ActiveCells, (self.size[0], self.size[1], self.size[2]), order="F")

    def read_prop(self, fname, prop, add=True, mult=1):
        """
        Reads input property from .dat file

        Parameters
        ----------
        fname: str
            Input file
        prop : str
            CMG keyword to read
        add : bool
            If grid geometry has been constructed, optionally add
            property data to grid as cell data with True
            False will return array_like of property data
        mult : float
            Multiplicative factor for scaling data
            Regrid leaves it to the user to determine when
            data needs to be scaled

        Returns
        -------
        data : array_like
            Property following natural cell ordering
        """
        print('Reading ' + prop + ' input')
        typeVal = None
        val = 0
        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == prop:
                        if len(item) >= 2:
                            if item[1] == "*CON":
                                val = float(item[2])
                                typeVal = '*CON'
                            elif item[1] == '*EQUALSI' or item[1] == 'EQUALSI':
                                attr_I = prop[:-1] + 'I'
                                # Change 'PERMJ' to be the keyword that identifies the end of attribute section
                                data = self.read_prop(fname, attr_I, add=False, mult=mult)
                                if len(item) == 4:
                                    op = item[2]
                                    if op == '*':
                                        data *= float(item[3])
                                    elif op == '/':
                                        data /= float(item[3])
                                    elif op == '+':
                                        data += float(item[3])
                                    elif op == '-':
                                        data -= float(item[3])
                            elif item[1] == 'ALL':
                                typeVal = 'ALL'
                        break

            if typeVal == 'ALL':
                data = []
                count = 0
                for line in fp:
                    item = line.split()
                    for attr in item:
                        if "*" in attr:
                            item = attr.split("*")
                            for i in range(0, int(item[0])):
                                data.append(float(item[1]))
                                count += 1
                        else:
                            data.append(float(attr))
                            count += 1
                    # If true, all values have been read
                    if count == self.size[0] * self.size[1] * self.size[2]:
                        data = np.array(data)
                        data = np.reshape(data, (self.size[0], self.size[1], self.size[2]), order="F")
                        break
            elif typeVal == '*CON':
                data = np.full((self.size[0], self.size[1], self.size[2]), val)

        if add:
            self.add_data(data, prop)
            self.out_props[prop] = data
        return data

    def read_ext_prop(self, fname, prop_title, mult=1):
        """
        Reads a property from an external file denoted by INCLUDE in .DAT
        Skips keywords at top of file, handles *MOD keyword at end of file

        Parameters
        ----------
        fname: str
            Input file
        prop_title : str
            Name for the data
        mult : float
            Multiplicative factor for scaling data
            Regrid leaves it to the user to determine when
            data needs to be scaled

        Returns
        -------
        data : array_like
            Property following natural cell ordering
        """
        print('Reading ' + prop_title + ' input')
        data = []
        count = 0
        modify = False
        with open(fname, "r") as fp:
            for line in fp:
                if not line[:1].isdigit():
                    if line.startswith('*MOD'):
                        modify = True
                    continue # it's a keyword
                item = line.split()
                if modify:
                    i = int(item[0])-1
                    j = int(item[1])-1
                    K = [int(x)for x in item[2].split(':')]
                    value = float(item[-1])
                    for k in range(K[0]-1,K[1]):
                        data[k,j,i] = value
                    break
                for attr in item:
                    if "*" in attr:
                        item = attr.split("*")
                        for i in range(0, int(item[0])):
                            data.append(float(item[1]) * mult)
                            count += 1
                    else:
                        data.append(float(attr) * mult)
                        count += 1
                # If true, all values have been read
                if count == self.size[0] * self.size[1] * self.size[2]:
                    data = np.array(data)
                    data = np.reshape(data, (self.size[2], self.size[1], self.size[0]), order="C")
                    continue
        self.add_data(data, prop_title)
        self.out_props[prop_title] = data

    def add_data(self, d, prop_title):
        """
        Adds data to a structured vtk grid

        Parameters
        ----------
        d : numpy_array
            Dimensions should match that of the grid adding to
        prop_title : str
            Name for the data
        """
        ac = vtk.vtkDoubleArray()
        ac.SetName(prop_title)
        for iac in d.flatten(order='F'):
            ac.InsertNextTuple1(iac)
        self.Grid.GetCellData().AddArray(ac)

    # Populates entire K-layer with val (for reading .out property)
    def _k_full(self, val):
        kLayer = dict((i, []) for i in range(self.size[1]))
        for j in range(self.size[1]):
            iRow = []
            for i in range(self.size[0]):
                iRow.append(val)
            kLayer[j] = iRow
        return kLayer

    # TODO: this is very similar to read_outputs, and is only used by refine
    def read_grid(self, data, dims):
        if dims[2] == 1:
            dims.pop()
        # layers = {}
        layers = np.zeros(dims)
        propIdxs = []
        I = None
        J = None
        K = '1'
        for line in data:
            item = line.split()
            if item[0] == 'All':
                # kKeys = np.arange(dims[2])
                # grid = dict((el, {}) for el in kKeys)
                # for k in range(dims[2]):
                #     kLayer = self._k_full(item[3])
                #     grid[k] = kLayer
                # self.out_props[attr_title][time] = grid
                # continue

                # numpy.full produces weird output when z axis is 1
                # TODO generalize this
                return np.full(dims, item[3])

            elif item[0] == 'Plane':
                K = item[3]
                if len(item) > 4:
                    if item[4] == 'All':
                        kLayer = self._k_full(item[7])
                        layers[int(K)-1] = kLayer
                else:
                    I = None
                    # layers[K] = {}
            if item[0] == 'I':
                # if K == '1' and I is None:
                    # layers[K] = {}
                J = None
                propIdxs = []
                I = item[2:]
                prevDigit = False
                for i in range(len(line)):
                    if line[i].isdigit():
                        if not prevDigit:
                            propIdxs.append(i)
                            prevDigit = True
                    else:
                        prevDigit = False
            # Check if there are any missing values in J line
            skipItem = []
            if item[0] == 'J=':
                JIdx = item[1]
                J = item[2:]
                # if JIdx not in layers[K].keys():
                    # layers[K-1][JIdx-1] = []
                for i in range(len(propIdxs)):
                    if line[propIdxs[i]] == ' ':
                        skipItem.append(i)
                numSkips = 0
                for i in range(len(I)):
                    if i in skipItem:
                        # layers[K][JIdx].append('NULL')
                        layers[int(JIdx)-1][i] = 'NULL'
                        numSkips += 1
                    else:
                        layers[int(JIdx)-1][i] = J[i - numSkips]
                        # layers[K][JIdx].append(J[i - numSkips])
            # Put entire grid worth of property in dictionary for current time step
            if I is not None and J is not None:
                # if int(I[-1]) == dims[0] and int(JIdx) == dims[1] and int(K) == dims[2]:
                if int(I[-1]) == dims[0] and int(JIdx) == dims[1]:
                    return layers

    def _is_float(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def _str_to_float(self, s):
        if not s[0].isdigit():
            s = -1
        else:
            if self._is_float(s):
                s = float(s)
            else:
                s = -2
        return s

    def read_outputs(self, fname, out_props, refined_blocks=None):
        """
        Reads per-cell output properties for every time step from .OUT file
        Properties available for reading determined by *OUTPRN *GRID
        If a cell property is empty in .OUT, then this will set it to null

        Parameters
        ----------
        fname : str
            File containing output properties (.OUT)
        out_props : array_like of array_like of str
            Specify multiple time-series properties to read from GEM .OUT file
            Sub-arrays should contain two strings:
                str 1 : An output property name as it appears at a timestep declaration
                        in a .OUT file (Time = 0 ... property name)
                str 2 : A custom name for that property- what it will be named if exported
                        to Numpy or VTK files
            ex)  [['Pressure (kpa)', 'Pressure'], ...]
        refined_blocks : dict
            Refinement blueprint from get_refined_blocks
            Only supports reading of 2D refinements
        """
        for prop in out_props:
            attr_name = prop[0].replace(" ", "").strip()
            attr_title = prop[1]
            self.out_props[attr_title] = {}
            if not hasattr(self, 'times'):
                self.times = []
            build = False
            found = False

            print('Reading ' + attr_title + ' output')
            with open(fname, "r") as fp:
                for line in fp:
                    item = line.split()
                    if len(item) > 0:
                        # Find current time step
                        if item[0] == 'Time':
                            time = item[2]
                            continue
                        attr = line.replace(" ", "").strip()

                        # Locate attribute name
                        if attr == attr_name:
                            build = True
                            found = True
                            layers = [[[] for j in range(self.size[1])] for i in range(self.size[2])]
                            I = None
                            J = None
                            K = '1'
                            continue

                        if build:
                            if item[0] == 'All':
                                v = self._str_to_float(item[3])
                                grid = np.full((self.size[2], self.size[1], self.size[0]), v, order='F')
                                self.out_props[attr_title][time] = grid
                                build = False
                                continue

                            if item[0] == 'Plane':
                                K = item[3]
                                if len(item) > 4:
                                    if item[4] == 'All':
                                        v = self._str_to_float(item[-1])
                                        k_layer = np.full((self.size[1], self.size[0]), v, order='F')
                                        layers[int(K)-1] = k_layer
                                        I = [self.size[0]]
                                        J = []
                                        JIdx = self.size[1]
                                else:
                                    I = None
                                continue

                            # Detects line indices under which 'I' properties are printed
                            if item[0] == 'I':
                                J = None
                                prop_i = []
                                I = item[2:]
                                prev_digit = False
                                for i in range(len(line)):
                                    if line[i].isdigit():
                                        if not prev_digit:
                                            prop_i.append([i])
                                        else:
                                            prop_i[-1].append(i)
                                        prev_digit = True
                                    else:
                                        prev_digit = False
                                continue

                            # Check if there are any missing values in J line
                            skip_i = []
                            if item[0] == 'J=':
                                JIdx = item[1]
                                J = item[2:]
                                for i in range(len(prop_i)):
                                    skip = True
                                    for j in prop_i[i]:
                                        if line[j] != ' ':
                                            skip = False
                                    if skip:
                                        skip_i.append(i)
                                n_skip = 0
                                for i in range(len(I)):
                                    if i in skip_i:
                                        # layers[int(K)-1][int(JIdx)-1].append('NULL')
                                        layers[int(K) - 1][int(JIdx) - 1].append(-1)
                                        n_skip += 1
                                    else:
                                        if J[i - n_skip][-1] == 'r':
                                            J[i - n_skip] = J[i - n_skip][:-1]

                                        v = self._str_to_float(J[i - n_skip])
                                        layers[int(K) - 1][int(JIdx) - 1].append(v)
                                        # layers[int(K)-1][int(JIdx)-1].append(J[i - n_skip])
                                continue

                            # Put entire grid worth of property in dictionary for current time step
                            if I is not None and J is not None:
                                if int(I[-1]) == self.size[0] and int(JIdx) == self.size[1] and int(K) == self.size[2]:
                                    if build:
                                        self.out_props[attr_title][time] = layers
                                        build = False

                                        # Add local grid refinements
                                        if refined_blocks:
                                            refined = self.refine_outputs(fp, copy.deepcopy(refined_blocks))
                                            for c in list(refined.keys()):
                                                idx = c.split(',')
                                                self.out_props[attr_title][time][idx[2]][idx[1]][int(idx[0])-1] = refined[c]
                if not found:
                    print(attr_title + ' was not found in output file!')
                    del self.out_props[attr_title]
                else:
                    if len(self.times) == 0:
                        self.times = list(self.out_props[attr_title].keys())

    # TODO: lightly tested, may have issues- only works with 2D
    def get_refined_blocks(self, fname):
        """
        To be used in grids with local grid refinement
        Inserts refined outputs into an LGR cell
        This reads the refinement blueprint for each cell, which can be fed to refine_outputs

        Parameters
        ----------
        fname : str
            Input file

        Returns
        -------
        refine_blocks : dict
            {cell idx: subdivision info, ...} for REFINE keywords in fname
        """
        refine_blocks = {}
        subgrid = []
        reading = False
        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if 'REFINE' in item[0]:
                        reading = True
                        if 'INTO' in line:
                            subgrid = [int(n) for n in item[-3:]]
                        refine_blocks[item[1]] = copy.deepcopy(subgrid)
                    else:
                        if reading:
                            break
        return refine_blocks

    # TODO: lightly tested, may have issues- only works with 2D
    def refine_outputs(self, fp, refined_blocks):
        """
            Refined grid in 39,27,81   Refined grid in 10,28,81
            ------------------------   --------------------------------
            All values are 0.113            I =  1        2        3
                                       J=  1   0.182    0.182    0.182
                                       J=  2   0.182    0.182    0.182
                                       J=  3   0.182    0.182    0.182

            Above is an example .OUT refinement
            We assume that the number of dashes is >= number of data items per line below
            This allows us to split data blocks based on number of dashes

            Currently, only 2D refiniements are supported

            Parameters
            ----------
            fp : textIOWrapper
                File object to read lines from
                ex) with open(fname, "r") as fp:
                        # move file pointer
                        ...
                        refine_outputs(fp, refined_blocks)
            refined_blocks : dict
                Refinement blueprint from get_refined_blocks
        """
        find_refine = True
        count_dash = False
        split_points = [0]
        block = []
        blocks_to_process = list(refined_blocks.keys())
        for line in fp:
            # Search for refined keywords
            if find_refine:
                idxs = line.strip('\n').split('Refined grid in')
                if len(idxs[0]) > 0:
                    idxs = [idx.strip(' ') for idx in idxs]
                    while '' in idxs:
                        idxs.remove('')
                    count_dash = True
                    find_refine = False
                    continue
            # count number of dashes to determine split points
            if count_dash:
                p = len(line) - len(line.lstrip(' '))
                dashes = line.strip('\n').split()
                for i in range(len(dashes) - 1):
                    p += len(dashes[i]) + 3
                    split_points.append(p)
                    # split_points.append(len(dashes[i]) + p)
                    # Assumes there are 3 spaces inbetween data blocks
                    # p += len(dashes[i]) + 3
                block = [[] for i in range(len(split_points))]
                count_dash = False
                continue
            # read data between split points
            if not find_refine and not count_dash:
                if len(line.strip('\n')) > 0:
                    # Use number of dashs to split lines of refinement data
                    data = [line[i:j] for i, j in zip(split_points, split_points[1:] + [None])]
                    for i,d in enumerate(data):
                        if len(d) > 0:
                            if 'I' in d:
                                d = d.rstrip()
                            else:
                                d = d.strip()
                            block[i].append(d)
                # refinement lines has been read, process data
                else:
                    for i,d in enumerate(block):
                        dims = refined_blocks[idxs[i]]
                        refined_blocks[idxs[i]] = self.read_grid(d, dims)
                        blocks_to_process.remove(idxs[i])
                    if len(blocks_to_process) == 0:
                        return refined_blocks
                    find_refine = True
                    split_points = [0]

    def get_wells(self, fname):
        """
        This should be ran before any other well operations are performed
        Builds a Wells object that contains dict of well init properties

        Parameters
        ----------
        fname : str
            File containing well init information

        Returns
        -------
        Wells : Wells object
            Well info read here assigned to Wells.wells
        """
        getLoc = False
        wells = []
        well = {'NAME': None, 'TYPE': None, 'OP_MODE': [], 'CON_TYPE': [], 'CON_VAL': [], 'LOC': None}
        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    keyword = item[0].strip('*')
                    if getLoc:
                        if item[0][0] == '*' and item[0][1] == '*':
                            continue
                        if ':' in item[2]:
                            K = [int(x)for x in item[2].split(':')]
                            well['LOC'] = (int(item[0]), int(item[1]), int(K[-1]))
                        else:
                            well['LOC'] = (int(item[0]), int(item[1]), int(item[2]))
                        getLoc = False
                    elif keyword == 'WELL':
                        # Add the previous well to the grid
                        if well['NAME'] is not None:
                            wells.append(copy.deepcopy(well))
                            well = {'NAME': None, 'TYPE': None, 'OP_MODE': [], 'CON_TYPE': [], 'CON_VAL': [],
                                    'LOC': None}
                        well['NAME'] = item[1].strip("'")
                    elif keyword == 'INJECTOR':
                        well['TYPE'] = 'INJ'
                    elif keyword == 'PRODUCER':
                        well['TYPE'] = 'PRO'
                    elif keyword == 'OPERATE':
                        well['OP_MODE'].append(item[1].strip('*'))
                        well['CON_TYPE'].append(item[2].strip('*'))
                        well['CON_VAL'].append(item[3].strip('*'))
                    elif keyword == 'PERF':
                        getLoc = True
            wells.append(well)
        return self.Wells(wells, self.Grid, self.times, self.out_dir)

    # TODO: kind of awkward, but has its use cases for visualizing wells
    def well_layers_vtk(self, wells):
        """
        NOTE: kind of awkward, but has its use cases for well viz
        Will probably be changed or removed at some point

        Adds well type/operational constraints to grid for vtk visualizations
        Labels injector well cells as 1, producer well cells as -1
        In Paraview, use Threshold filter to see these cells

        Parameters
        ----------
        wells : array_like of dict
            After Grid.get_wells(), Wells.wells
        """
        all_wells = np.zeros(self.size[0] * self.size[1] * self.size[2])
        well_data_layers = {'ALL': all_wells, 'INJ': {}, 'PRO': {}}
        for well in wells:
            i = 0
            loc = well['LOC']
            idx = ((self.size[0] * self.size[1]) * (loc[2] - 1)) + (self.size[0] * (loc[1] - 1)) + (loc[0] - 1)
            if well['TYPE'] == 'INJ':
                for con_type in well['CON_TYPE']:
                    if con_type not in well_data_layers['INJ'].keys():
                        well_data_layers['INJ'][con_type] = np.zeros(self.size[0] * self.size[1] * self.size[2])
                    else:
                        well_data_layers['INJ'][con_type][idx] = well['CON_VAL'][i]
                    i += 1
                well_data_layers['ALL'][idx] = 1
            else:
                for con_type in well['CON_TYPE']:
                    if con_type not in well_data_layers['PRO'].keys():
                        well_data_layers['PRO'][con_type] = np.zeros(self.size[0] * self.size[1] * self.size[2])
                    i += 1
                well_data_layers['ALL'][idx] = -1

        # Add data layers to grid
        data_layer = well_data_layers['ALL']
        self.add_data(data_layer, 'WELLS ALL')
        for con_type in well_data_layers['INJ'].keys():
            data_layer = well_data_layers['INJ'][con_type]
            self.add_data(data_layer, 'WELLS INJ ' + con_type)
        for con_type in well_data_layers['PRO'].keys():
            data_layer = well_data_layers['PRO'][con_type]
            self.add_data(data_layer, 'WELLS PRO ' + con_type)

    # TODO: move outside, generalize across simulators?
    class Wells:
        """
        Contains useful methods for working with wells
        Should not be called explicitely, use Grid.get_wells() to create
        """
        def __init__(self, wells, grid, times, out_dir):
            self.wells = wells
            self.grid = grid
            self.times = times
            self.out_dir = out_dir

        # Return list of well names
        def names(self):
            return [well['NAME'] for well in self.wells]

        # Get a single well by name
        def well(self, name):
            for well in self.wells:
                if well['NAME'] == name:
                    return well
            raise KeyError('Well name not recognized')

        # Compute center of cell given corner point coords
        def centroid(self, coords):
            return np.mean(coords, axis=0)

        # Get xyz coords of corner points defining cell
        # TODO: Duplicated from grid object
        def cell_verts(self, idx):
            coords = []
            xyz = [0, 0, 0]
            cell = self.grid.GetCell(idx[0]-1, idx[1]-1, idx[2]-1)
            p_ids = cell.GetPointIds()
            n_ids = p_ids.GetNumberOfIds()
            for n in range(n_ids):
                p = p_ids.GetId(n)
                self.grid.GetPoint(p, xyz)
                coords.append(copy.deepcopy(xyz))
            return coords

        def _out_order(self, fname):
            """
            Helper fn for readOutput
            In GEMFIELDSUMMARY blocks, not all wells may be listed
            This finds which ones are and returns an ordered list of their names
            """
            # t = 1
            orderDict = {}
            order = []
            readWells = False
            lastBlock = False
            addOrder = False
            with open(fname, "r") as fp:
                for line in fp:
                    item = line.split()
                    if readWells:
                        if lastBlock:
                            line = line.split('++')[0]
                            addOrder = True
                            lastBlock = False
                        item = list(map(str.strip, line.split('+')))
                        item = [e.split() for e in list(filter(None, item))]
                        order.extend([w[1] for w in item])
                        readWells = False
                        if addOrder:
                            orderDict[t] = order
                            order = []
                            addOrder = False
                            # t += 1
                    elif len(item) > 0:
                        head = ''.join(item[2:])
                        if 'GEMFIELDSUMMARY' in head:
                            t = item[1]

                        elif 'No.' in line and 'Name' in line and '+' in line:
                            if '++' in line:
                                lastBlock = True
                            readWells = True
                            next(fp)
                            continue
            return orderDict

        def build_cylinders(self, vtk_fname, zscale=1, radius=20):
            """
            Prepares Paraview script for automatically creating cylinders to represent wells
            Generated Paraview script can be executed in Paraview Python Shell after VTK grid
            has been loaded
            Grid.get_wells() must be called first to gather well info and create wells object
            Does not currently support wells with multiple perforations

            Parameters
            ----------
            vtk_fname : str
                Should be the same as that provided to export_grid(vtk_fname)
            zscale : float
                Should match any Z-scaling (planning to be) applied to the grid in Paraview
            radius : float
                Cylider radius

            Returns
            -------
            py_file : .py file
                Python/Paraview script for creating cylinders at well locations in Paraview
            """
            centers = []
            heights = []
            for w in self.wells:
                b = w['LOC']
                t = [b[0], b[1], 1]
                b_surf = self.cell_verts(b)
                t_surf = self.cell_verts(t)
                minZ = t_surf[0][2]
                b_v = self.centroid(b_surf)
                for i in range(len(b_surf)):
                    t_v = t_surf[i]
                    if t_v[2] < minZ:
                        minZ = t_v[2]
                cyl_height = abs(b_v[2] - minZ)
                cyl_center = [b_v[0], b_v[1], minZ + cyl_height/2]
                centers.append(cyl_center)
                heights.append(cyl_height)
            # Write center, height data to string
            c_s = '['
            for c in centers:
                c_s += '[' + ', '.join(map(str, c)) + '],\n'
            c_s = c_s[:-2] + ']\n\n'
            h_s = '['
            for h in heights:
                h_s += str(h) + ',\n'
            h_s = h_s[:-2] + ']\n\n'
            # Write Paraview script
            with io.open('well_cylinders_test.py', 'w', newline='\r\n') as f:
                f.write("from paraview.simple import *\n")
                f.write("paraview.simple._DisableFirstRenderCameraReset()\n\n")
                f.write("# To use, open grid vtk files in Paraview, then run this script in Python Shell\n")
                f.write("c = " + c_s)
                f.write("h = " + h_s)
                f.write("for i in range(len(h)):\n")
                f.write("    cylinder = Cylinder()\n")
                f.write("    fn = FindSource('" + vtk_fname + "_*')\n")
                f.write("    rv = GetActiveViewOrCreate('RenderView')\n")
                f.write("    tf = Transform(Input=cylinder)\n")
                f.write("    tf.Transform = 'Transform'\n")
                f.write("    tf.Transform.Translate = [c[i][0], c[i][1], c[i][2] * " + str(zscale) + "]\n")
                f.write("    tf.Transform.Rotate = [90.0, 0.0, 0.0]\n")
                f.write("    tf.Transform.Scale = [1.0, " + str(zscale) + ", 1.0]\n")
                f.write("    transformDisplay = Show(tf, rv, 'GeometryRepresentation')\n")
                f.write("    cylinder.Radius = " + str(radius) + "\n")
                f.write("    cylinder.Height = h[i]\n")
                f.write("    rv.Update()")

        def read_output(self, fname, keys, subkeys):
            """
            Reads well response/output information for all time steps
            The properties being read are located at the 'GEM FIELD SUMMARY'
            at each time step

            Parameters
            ----------
            fname : str
                File containing well output properties (.OUT)
            keys : array_like of str
                These are the property titles, which act as a header for the
                properties (subkeys) which are to be read
                These can be gathered by locating a 'GEM FIELD SUMMARY' section
                ex) ['Well Pressures', 'Inst Surface Production Rates']
            subkeys : array_like of array_like of str
                These are the actual properties for which to read data
                These can be gathered by locating a 'GEM FIELD SUMMARY' section
                A subarray index should match its corresponding key index
                    (make sure this order is correct!)
                ex) This example goes along with the keys example above:
                    [['Bottom Hole', 'Drawdown'], ['Oil', 'Water', 'Gas']]

            Returns
            -------
            well_out : dict
                Output data for all keys/subkeys for all wells for all time steps
            """
            print('Reading well outputs')
            k = 0
            sk = 0
            # tID = 0
            build_keys = False
            build_subkeys = False
            cur_key = None
            cur_block = 1
            n_blocks = 0
            order = self._out_order(fname)
            # Initialize output dictionary
            well_out = {}
            for i,key in enumerate(keys):
                well_out[key] = {}
                for subkey in subkeys[i]:
                    well_out[key][subkey] = {}
                    # for t in range(len(self.times) - 1):
                    #     well_out[key][subkey][t+1] = []
            with open(fname, "r") as fp:
                for line in fp:
                    item = line.split()
                    if len(item) > 0:
                        # Find current time step
                        if item[0] == 'TIME:':
                            head = ''.join(item[2:])
                            if 'GEMFIELDSUMMARY' in head:
                                build_keys = True
                                # tID += 1
                                tID = item[1]
                                n_wells = len(order[tID])
                                n_blocks = math.ceil(n_wells / 4)
                                continue
                        # Assume that keywords are ordered as they appear in .out
                        if build_keys:
                            # Current block has been read, move to next one
                            if k == len(keys):
                                cur_block += 1
                                k = 0
                            if keys[k] in line:
                                cur_key = keys[k]
                                build_keys = False
                                build_subkeys = True
                                continue
                        elif build_subkeys:
                            if sk == len(subkeys[k]):
                                build_subkeys = False
                                build_keys = True
                                k += 1
                                sk = 0
                                continue
                            for subkey in subkeys[k]:
                                if subkey in line:
                                    if cur_block == n_blocks:
                                        line = line.split('++')[0]
                                    item = list(map(str.strip, line.split('+')))
                                    item = list(filter(None, item))
                                    if tID not in well_out[cur_key][subkey]:
                                        well_out[cur_key][subkey][tID] = item[1:]
                                    else:
                                        well_out[cur_key][subkey][tID].extend(item[1:])
                                    sk += 1
            # Attach well names to well outputs
            for key in well_out:
                for subkey in well_out[key]:
                    # for t in range(1, len(self.times)):
                    for t in well_out[key][subkey]:
                        well_out[key][subkey][t] = {k:v for k,v in zip(order[t], well_out[key][subkey][t])}
            return well_out

    def _prep_vtk(self, data, prop, propIds):
        vtk_data = data.flatten(order='C')
        ac = vtk.vtkDoubleArray()
        ac.SetName(prop)
        for iac in vtk_data:
            ac.InsertNextTuple1(iac)
        id = self.Grid.GetCellData().AddArray(ac)
        propIds.append(id)
        return propIds

    def export_grid(self, vtk_fname='GRID', toVTK=True, toNumpy=True):
        """
        Export grid information to Numpy or VTK files after all data has been read

        Parameters
        ----------
        vtk_fname : str
            Output file name for generated VTK files
        toVTK : bool
            Optionally export grid data to VTK files
        toNumpy : bool
            Optionally export grid data to Numpy files
        """
        print('Exporting grids')
        tID = 0
        # Start by exporting input properties (from read_prop() or read_ext_prop())
        # In VTK files, these props will only be visible at only the first timestep
        dp = []
        propIds = []
        for prop in self.out_props:
            if type(self.out_props[prop]) is not dict:
                data = np.array(self.out_props[prop])
                # Save to Numpy
                if toNumpy:
                    self.export_prop(data, prop, tID)
                # Add property data to vts structured grid
                if toVTK:
                    propIds = self._prep_vtk(data, prop, propIds)
                    self._check_out('vtk')
            else:
                dp.append(prop)

        # Export time-series output properties (from read_out_props())
        for t in self.times:
            for prop in self.out_props:
                if prop in dp:
                    data = np.array(self.out_props[prop][t], order='F')
                    # Save to Numpy
                    if toNumpy:
                        self.export_prop(data, prop, tID)
                    # Add property data to vts structured grid
                    if toVTK:
                        propIds = self._prep_vtk(data, prop, propIds)
            # Save to VTK
            if toVTK:
                if tID == 0:
                    self._check_out('vtk')
                self.exportVTK(os.path.join(self.out_dir, 'vtk', vtk_fname + str(tID)))
                for id in propIds:
                    self.Grid.GetCellData().RemoveArray(id)
            tID += 1
            propIds = []

    def export_prop(self, d, title, t):
        """
        Exports a Numpy array of a grid property for a given timestep

        d : array_like
            Grid data to be saved
        title : str
            Name for the data
        t : int
            Current timestep id

        """
        self._check_out(title)
        np.savez_compressed(os.path.join(self.out_dir, title, title + '_' + str(t)), d)

    def export_wells(self, w, title):
        """
        Exports well dictionaries as a Numpy file

        Parameters
        ----------
        w : dict
            Well data can be fetched from two places:
                Wells.wells for well location, type, constraints
                Wells.read_outputs() for well response/outputs over time
        title : str
            Name for the data
        """
        self._check_out(title)
        np.savez_compressed(os.path.join(self.out_dir, title, title), w)

# TODO: under development
class DynamicGrid:
    """
    UNDER DEVELOPMENT
    """
    def __init__(self, wells, grid, times, out_dir):
        self.wells = wells
        self.grid = grid
        self.times = times
        self.out_dir = out_dir

    def evolve(self, prop, cell_filter_fn, prop_change_fn):
        pass