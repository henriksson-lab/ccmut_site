#!/usr/bin/env python3
import sqlite3
import re
import pandas as pd
import numpy as np


def decompressMSARLE(s):
	return ''.join(['-'*int(x) if x.isnumeric() else x
		for x in re.split(r'(\d+)', s)])




listedit = pd.read_csv("listedit.csv")

numedit = listedit.shape[0]

conn = sqlite3.connect("file:alignments.sqlite?mode=ro", uri=True)
for i in range(0,numedit):
	sqldata = pd.read_sql_query("SELECT strain, wbid, alignment from alignment_cropped where strain=? and wbid=?", conn,params=[listedit["ccid"][i], listedit["wbid"][i]])

	al=decompressMSARLE(sqldata["alignment"].tolist()[0])
	print()
	print()
	print()
	print(al)
	print()
	print(listedit["ccid"][i] + "     " + listedit["wbid"][i])
	print()

conn.close()
