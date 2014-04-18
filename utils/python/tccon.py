#
# Python library to access TCCON data
#
# Dietrich Feist, Max-Planck-Institut for Biogeochemistry, Jena, Germany
# 2013-10-16
#

#
# General class to read TCCON data files
#
class tccon_file:
    """General class to read TCCON data files"""
    
    # Separator string (use None for whitespace)
    sep = None
    
    # Initialize class object
    def __init__(self, file_name):
        # Initialize parameters
        self.file_name = file_name
        
        # Read whole file
        with open(self.file_name, 'r') as self.file:
            self.line1_par = self.read_line1()
            self.header = self.read_header()
            self.init_columns_units()
            self.data = self.read_data()
        
    # Read 1st line that contains several integers
    def read_line1(self):
            # Read first line
            line1_par = self.file.readline().split()
            
            # Try to convert all values to integers
            for i in range(len(line1_par)):
                try:
                    line1_par[i] = int(line1_par[i])
                except ValueError:
                    line1_par = []
                    break 
            return line1_par
        
    # Read file header -> list of one string per line
    def read_header(self):
        header = []
        if self.line1_par:
            # Go to beginning of header
            self.file.seek(0)
            self.file.next()                
            # Append each line to header
            for i in range(self.line1_par[0] - 1):
                header.append(self.file.next().rstrip())
        return header

    # Init column and unit names (default)
    def init_columns_units(self):
        self.columns = []
        for i in range(self.line1_par[1]):
            self.columns.append('%d' % i)
        self.units = self.line1_par[1] * ['']
                
    # Read data from file
    def read_data(self):        
        if self.line1_par:
            # Skip to beginning of data section
            self.file.seek(0)
            for i in range(self.line1_par[0]):
                self.file.next()
                
            # Initialize data structure
            data = {}
            for column in self.columns:
                data[column] = []
                
            # Convert each data line
            l = 0
            for line in self.file:                
                row = self.par_split(line)
                l += 1
                
                # Check for format problems
                if len(row) != len(self.columns):
                    print "Format error in data block line #%d: expected %d columns, received %d columns!" % (l, len(self.columns), len(row))
                    continue
                
                # Append row data
                for i in range(len(row)):
                    # Try to convert value to int if possible
                    if row[i].strip().isdigit():
                        val = int(row[i])
                    else:
                        # Try to convert value to float or leave as string
                        try:
                            val = float(row[i])
                        except ValueError:
                            val = row[i]
                    # Append value      
                    data[self.columns[i]].append(val)
                    
            # Convert columns to numpy arrays
            import numpy
            for col in data.keys():
                data[col] = numpy.array(data[col])
                
        return data
    
    # Split line into list of parameters
    def par_split(self, line):
        # line:  string read from a data file
        par = line.strip().split(self.sep)
        return par

#
# Special class to read oof (official output) files
# Column names provided in last header line, no units
#
class oof_file(tccon_file):
    """Special class to read TCCON oof (official output) files"""
    
    # Init column and unit names to values provided in header
    def init_columns_units(self):
        self.columns = self.par_split(self.header[-1])
        self.units = self.line1_par[1] * ['']

#
# Special class to read col files (similiar to oof)
#
class col_file(oof_file):
    """Special class to read col files"""

#
# Special class to read eof files (similiar to oof)
#
class eof_file(oof_file):
    """Special class to read eof files"""

#
# Special class to read gop (sunrun) files (similiar to oof)
#
class gop_file(oof_file):
    """Special class to read gop (sunrun) files"""

#
# Special class to read grl (runlog) files (similiar to oof)
#
class grl_file(oof_file):
    """Special class to read grl (runlog) files"""

#
# Special class to read ray files (similiar to oof)
#
class ray_file(oof_file):
    """Special class to read ray files"""

#
# Special class to read tav files (similiar to oof)
#
class tav_file(oof_file):
    """Special class to read tav files"""

#
# Special class to read tsw files (similiar to oof)
#
class tsw_file(oof_file):
    """Special class to read tsw files"""

#
# Special class to read vav files (similiar to oof)
#
class vav_file(oof_file):
    """Special class to read vav files"""

#
# Special class to read vsw files (similiar to oof)
#
class vsw_file(oof_file):
    """Special class to read vsw files"""

#
# Special class to read map (a priori) files
#
class map_file(tccon_file):
    """Special class to read map (a priori) files"""
    
    # map files are comma-separated
    sep = ','
    
    # Init column and unit names to values provided in header
    def init_columns_units(self):
        self.units = self.par_split(self.header[-1])
        self.columns = self.par_split(self.header[-2])
