#!/usr/bin/env python
# using Python 2

import sys
import os

usage="""
Usage:  txt2html.py predzinc-resultfile -o out-html-file
Options:
  -h|--help       : print this help message and exit

Created 2021-03-09, updated 2021-03-09, Nanjiang Shu

"""

def PrintHelp():# {{{
    print usage# }}}
def CreateNiceHTMLFile(predfile, nicehtmlfile):# {{{
# generate HTML file from the prediction result
    try:
        fpout = open(nicehtmlfile, "w")
    except IOError:
        fpout = None

    if fpout:
        fpout.write("<table class=\"frag1dresult\">\n")
        fpout.write("<tr>\n")
        fpout.write("<th>No.</th>\n")
        fpout.write("<th>AA</th>\n")
        fpout.write("<th>Sec</th>\n")
        fpout.write("<th>ConfSec</th>\n")
        fpout.write("<th>S8</th>\n")
        fpout.write("<th>ConfS8</th>\n")
        fpout.write("<th>S3</th>\n")
        fpout.write("<th>ConfS3</th>\n")
        fpout.write("</tr>\n")
    lines = open(predfile, "r").read().split("\n")
    for line in lines:
        if not line or line[0] == "#" or line[0] == "/":
            continue
        strs = line.split()
        if len(strs) == 8:
            sec = strs[2]
            s3 = strs[6]
            if sec == "H":
                bgcolor_sec = "red"
            elif sec == "S":
                bgcolor_sec = "yellow"
            else:
                bgcolor_sec = ""

            if s3 == "H":
                bgcolor_s3 = "red"
            elif s3 == "S":
                bgcolor_s3 = "yellow"
            elif s3 == "T":
                bgcolor_s3 = "lightblue"
            else:
                bgcolor_s3 = ""

            if fpout:
                fpout.write("<tr>\n")
                fpout.write("<td align=\"center\">%s</td>\n"%(strs[0])) #No
                fpout.write("<td align=\"center\">%s</td>\n"%(strs[1])) #AA
                fpout.write("<td bgcolor=\"%s\" align=\"center\">%s</td>\n"%(bgcolor_sec, strs[2])) #Sec
                fpout.write("<td align=\"center\">%s</td>\n"%(strs[3])) #ConfSec
                fpout.write("<td align=\"center\">%s</td>\n"%(strs[4])) #S8
                fpout.write("<td align=\"center\">%s</td>\n"%(strs[5]))  #ConfS8
                fpout.write("<td bgcolor=\"%s\" align=\"center\">%s</td>\n"%(bgcolor_s3, strs[6])) #S3
                fpout.write("<td align=\"center\">%s</td>\n"%(strs[7]))  #ConfS3
                fpout.write("</tr>\n")
    if fpout:
        fpout.write("</table>\n")
    if fpout:
        fpout.close()# }}}
def main():#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    predfile = ""
    out_htmlfile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            predfile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                out_htmlfile = sys.argv[i+1]
                i += 2
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            predfile = argv[i]
            i += 1


    if predfile != "" and out_htmlfile != "":
        CreateNiceHTMLFile(predfile, out_htmlfile)
        return 0
    else:
        print >> sys.stderr, "Either predfile or out_htmlfile not set. Exit."
        PrintHelp()
        return 1
#}}}
if __name__ == '__main__' :# {{{
    sys.exit(main())# }}}
