#Python function to convert Horn Radio Astronomy files from and to .csv files
#HISTORY
#21Oct13 GIL initial version based on rasnames.py
import os
import glob
import rasnames

def toCsv( names, doDebug=False):
    """
    toCsv convert .ast and .hot files to .csv files
    Inputs
    names - list of directories and files names
    doDebug indicates showing debugging info
    The function only searches one directory down from the current directory.
    """
    typea = str(".ast")   
    typeb = str(".hot")
    # based on older code flag finding both types of observations
    
    fullnames, nfull = rasnames.splitNames( names, typea, typeb, doDebug=doDebug)
    
    if doDebug:
        print("Number of full names: %d" % (nfull))
        print("Full Names:")
        print(fullnames)

    outcount = 0
    for aname in fullnames:
#        outfile.write('# File: ' + outname + '\n')
        try:
            if os.path.isfile(aname):
                fin = open(aname, 'r')
            else:
                print("Not a valid file name: %s" % (aname))
                continue
        except:
            print("Can Not open File: %s" % (aname))
            continue
        
# read the whole file into a single variable, which is a list of every row of the file.
        inlines = fin.readlines()
        fin.close()

        # create the output name from the input name
        outparts = aname.split(".")
        # replace the extension
        if len(outparts) > 1:
            outname = outparts[0] + ".csv"
        else:
            outname = outparts + ".csv"
            print("Unusual input file name: %s" % (aname))
            print("Outname: %s" % (outname))
        fout = open(outname, "w")
        # now prepare to return the new names
        if outcount == 0:
            outnames = [outname]
        else:
            outnames.append(outname)
            
        outcount = outcount + 1
# now parse the lines one at a time and produce a CSV version file
        for aline in inlines:
            # get rid of extra whitespace
            aline = aline.strip()
            # if this a keyword or comment
            if aline[0] == "#":
                # remove the first character
                aline = aline[1:]
                lineparts = aline.split("=")
                # if a keyword = value line
                if len(lineparts) > 1:
                    keyword = lineparts[0]
                    keyword = keyword.strip()
                    value = lineparts[1]
                    value = value.strip()
                    outline = "#,%-8s, %s" % (keyword, value)
                    # if a # in the line text
                    if len(lineparts) > 2:
                        outline = outline + "=" + lineparts[2]
                        if len(lineparts) > 3:
                            outline = outline + "=" + lineparts[3]                            
                else:
                    outline = "#,        , " + aline.strip()
            else:  # else no beginning #, must be data
                lineparts = aline.split(" ")
                nline = len(lineparts)
                if nline < 2:
                    print("Unusual line with only one value: %s" % (aline))
                    outline = aline
                else:
                    outline = lineparts[0]
                    for iii in range(nline-1):
                        outline = outline + ",  " + lineparts[iii+1]
                # end else data

            outline = outline + "\n"
            fout.write(outline)
            # end for all lines 
        fout.close()
    # end for all files
    # end of toCsv
    return outnames, outcount
    
def toAst( names, doDebug=False):
    """
    toAst converts .csv files to .ast or .hot files
    Inputs
    names - list of directories and files names
    doDebug indicates showing debugging info
    The function only searches one directory down from the current directory.
    """
    typea = str(".csv")   
    typeb = str("")
    # based on older code flag finding both types of observations
 
    # transform directory names into complete file names for any matching type
    csvnames, nfull = rasnames.splitNames( names, typea, typeb, doDebug=doDebug)
    
    outcount = 0
    
    for aname in csvnames:
     
#        outfile.write('# File: ' + outname + '\n')
        try:
            if os.path.isfile(aname):
                fin = open(aname, 'r')
            else:
                print("Not a valid file name: %s" % (aname))
                continue
        except:
            print("Can Not open File: %s" % (aname))
            continue
        
# read the whole file into a single variable, which is a list of every row of the file.
        inlines = fin.readlines()
        fin.close()

        # create the output name from the input name
        outparts = aname.split(".")
        # replace the extension
        if len(outparts) > 1:
            outname = outparts[0] + ".ast"
        else:
            outname = outparts + ".ast"
            print("Unusual input file name: %s" % (aname))
            print("Outname: %s" % (outname))
        fout = open(outname, "w")
        elevation = 90.  # assume a .ast file, may discover this is a .hot file
        # now prepare to return the new names
# now parse the lines one at a time and produce a AST version file
        for aline in inlines:
            aline = aline.strip()
            # if this a keyword or comment
            if aline[0] == "#":
                # remove the first character
                aline = aline[1:]
                aline = aline.strip()
                # now if a leading comma, remove it.
                if aline[0] == ',':
                    aline = aline[1:]
                    aline = aline.strip()
                lineparts = aline.split(",")
                # if a keyword = value line
                if len(lineparts) > 1:
                    keyword = lineparts[0]
                    keyword = keyword.strip()
                    value = lineparts[1]
                    value = value.strip()
                    if keyword == "EL":
                        elevation = float(value)
                    # if no keyword, then just a comment
                    if keyword == "":
                        outline = "# %s" % (value)
                    else:  # else there is a keyword
                        outline = "# %-10s= %s" % (keyword, value)   
                    # if a # in the line text
                    if len(lineparts) > 2:
                        outline = outline + " # " + lineparts[2]
                else:
                    outline = "# " + aline.strip()
            else:  # else no beginning #, must be data
                lineparts = aline.split(",")
                nline = len(lineparts)
                if nline < 2:
                    print("Unusual line with only one value: %s" % (aline))
                    outline = aline
                else:
                    outline = lineparts[0].strip()
                    for iii in range(nline-1):
                        outline = outline + " " + lineparts[iii+1].strip()
                # end else data

            outline = outline + "\n"
            fout.write(outline)
            # end for all lines 
        fout.close()
        
        # if a hot observation, must move from .ast to .hot
        if elevation < 0.:
            hotname = outparts[0] + ".hot"
            os.rename(outname, hotname)
            # prepare to report the hotfile name
            outname = hotname
        # if starting a list of files, then init else append the name to list        
        if outcount == 0:
            outnames = [outname]
        else:
            outnames.append(outname)
            
        outcount = outcount + 1


    # end for all files
    # end of toAst()
    return outnames, outcount
