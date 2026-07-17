import numpy as np


class InterpolatedEstuaryVelocity:
    """Functor to interpolate from an unstructured list of inflow velocity data to AMR-Wind boundary points."""
    # Interpolate velocity component from the FVCOM data to the AMR-Wind grid.  We also pass the surface elevation and water
    # height fields.  Set the velocity to subterrain value if below the Earth surface or to the air value if above the water 
    # level.
    
    def __init__(self,interpolatorField,interpolatorElev,interpolatorWaterLev,subTerrainVal,airVal):
        self.interpolatorField = interpolatorField
        self.interpolatorElev = interpolatorElev
        self.interpolatorWaterLev = interpolatorWaterLev
        self.subTerrainVal = subTerrainVal
        self.airVal = airVal

    
    def __call__(self, xg, yg, zg, time):        
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape

        # Interpolate the field variable, the surface elevation, and water level to the boundary points
        # using the scipy interpolator passed here.
        field = self.interpolatorField(xg,yg,zg)
        surfElev = self.interpolatorElev(xg,yg)
        waterLev = self.interpolatorWaterLev(xg,yg)

        # Check to see if the point actually lies within the water.  If so, use the interpolated from
        # FVCOM data value; if below the surface, use the subterrain value; if in the air, use the air value.
        # Because the boundary points are i, j, k ordered, use the .flat to iterate through the flattened
        # array without actually flattening.
        field_ = field.flat
        zg_ = zg.flat
        surfElev_ = surfElev.flat
        waterLev_ = waterLev.flat
        
        for i in range(len(field_)):
            if (zg_[i] > waterLev_[i]):
                field_[i] = self.airVal
            elif (zg_[i] < surfElev_[i]):
                field_[i]= self.subTerrainVal


       #self.contPlot(yg,zg,field)

        return field

    
    #def contPlot(self,y,z,f):
       #print("Placeholder to plot contours")




        


class InterpolatedEstuaryDensityVOF:
    """Functor to interpolate from an unstructured list of inflow velocity data to AMR-Wind boundary points."""
    # Set the density or VOF field on the AMR-Wind grid.  We also pass the water height fields.  Set to air value if above
    # the water level or water level if below, even if inside terrain.
    
    def __init__(self,interpolatorWaterLev,waterVal,airVal):
        self.interpolatorWaterLev = interpolatorWaterLev
        self.waterVal = waterVal
        self.airVal = airVal

    
    def __call__(self, xg, yg, zg, time):        
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape

        # Interpolate the water level to the boundary points using the scipy interpolator passed here.
        waterLev = self.interpolatorWaterLev(xg,yg)

        # Check to see if the point actually lies within the water or air and set the value accordingly.
        field = self.airVal*np.ones(xg.shape)
        
        field_ = field.flat
        zg_ = zg.flat
        waterLev_ = waterLev.flat
        for i in range(len(field_)):
            if (zg_[i] < waterLev_[i]):
                field_[i] = self.waterVal


       #self.contPlot(yg,zg,field)

        return field

    
    #def contPlot(self,y,z,f):
       #print("Placeholder to plot contours")




        


class Interpolated2D:
    """Functor to interpolate from an unstructured list of inflow data.  Assumes data in a 2D Cartesian plane."""
    # Initialize the interpolated functor class.  Pass the scipy.interpolate
    # interpolator along with the coordinates of the bottom of the terrain and
    # the orientation.  The orientation corresponds to xlo = 0, ylo = 1, 
    # zlo = 2, xhi = 3, yhi = 4, and zhi = 5
    def __init__(self,interpolator,orientation,xBottom,yBottom,zBottom,subTerrainValue):
        self.interpolator = interpolator
        self.orientation = orientation
        self.xBottom = xBottom
        self.yBottom = yBottom
        self.zBottom = zBottom
        self.subTerrainValue = subTerrainValue

    
    def __call__(self, xg, yg, zg, time):        
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape

        if   ((self.orientation == 0) or (self.orientation == 3)):
            x1 = yg
            x2 = zg
            xb1 = self.yBottom
            xb2 = self.zBottom
        elif ((self.orientation == 1) or (self.orientation == 4)):
            x1 = xg
            x2 = zg
            xb1 = self.xBottom
            xb2 = self.zBottom
        elif ((self.orientation == 2) or (self.orientation == 5)):
            x1 = xg
            x2 = yg
            
        
        field = self.interpolator(x1,x2)

        if not ((self.orientation == 2) or (self.orientation == 5)):
            field_ = field.flat
            x_ = x1.flat
            z_ = x2.flat
            for i in range(len(field_)):
                zTerrain = np.interp(x_[i],xb1,xb2)
                if (z_[i] < zTerrain):
                    field_[i] = self.subTerrainValue

       #self.contPlot(yg,zg,field)

        return field

    
    #def contPlot(self,y,z,f):
       #print("Placeholder to plot contours")




        


class Constant:
    """Constant functor."""

    def __init__(self, constant):
        self.constant = constant

    def __call__(self, xg, yg, zg, time):
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape
        
        return self.constant * np.ones(xg.shape)







class Tabular1D:
    """Linearly interpolate from table functor."""

    def __init__(self, dir, table):
        self.dir = dir
        self.table = table
        self.coord = []
        self.field = []


    def __call__(self, xg, yg, zg, time):
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape
        
        coord = np.zeros(np.shape(xg))
        if ((self.dir == 0) or (self.dir == 'x')):
            self.coord = xg
        elif ((self.dir == 1) or (self.dir == 'y')):
            self.coord = yg
        elif ((self.dir == 2) or (self.dir == 'z')):
            self.coord = zg

            
        coordMin = np.min(self.coord)
        coordMax = np.max(self.coord)
        print(len(self.coord),coordMin,coordMax)

        # = np.zeros(np.shape(self.coord))

        self.field = np.interp(self.coord,self.table[0],self.table[1])
        
        return self.field

    def getCoord(self):
        return self.coord

    def getField(self):
        return self.field

    def getUnique(self):
        coordUnique, indices = np.unique(self.coord.flat, return_index=True)
        fieldUnique = self.field.flat[indices]
        return coordUnique, fieldUnique





class PoiseuilleFlow:
    """One dimensional parabolic Poiseuille functor."""

    def __init__(self, dir, D, U, profileFile):
        self.dir = dir 
        self.D = D
        self.U = U
        self.profileFile = profileFile

    
    def __call__(self, xg, yg, zg, time):        
        assert xg.shape == yg.shape
        assert xg.shape == zg.shape
        
        coord = np.zeros(np.shape(xg))
        if ((self.dir == 0) or (self.dir == 'x')):
            coord = xg
        elif ((self.dir == 1) or (self.dir == 'y')):
            coord = yg
        elif ((self.dir == 2) or (self.dir == 'z')):
            coord = zg
        
        field = np.zeros(np.shape(coord))
        
        field_ = field.flat
        coord_ = coord.flat
        for i in range(len(field_)):
            if ((coord_[i] >= 0.0) and (coord_[i] <= self.D)):
                field_[i] = 4.0 * self.U * ((coord_[i]/self.D) - (coord_[i]/self.D)**2)
            else:
                field_[i] = 0.0

        self.profile(coord, self.profileFile)

        return field

    
    def profile(self, coord, profileFile):


        coord = np.sort(np.unique(coord))
        field = np.zeros(np.shape(coord))

        fid = open(profileFile,'w')
        
        for i in range(len(coord)):
            if ((coord[i] >= 0.0) and (coord[i] <= self.D)):
                field[i] = 4.0 * self.U * ((coord[i]/self.D) - (coord[i]/self.D)**2)
            else:
                field[i] = 0.0

            fid.write(f"{coord[i]:.16e} {field[i]:.16e} 0.0 0.0 0.0\n")
        
        fid.close()