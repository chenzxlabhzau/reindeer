import MySQLdb
db = MySQLdb.connect("localhost","nazhang","nazhang_123","orthomcl_nazhang")
cursor = db.cursor()

f = open("/home/mwshi/project/xunlu/mcl/similarSequences.txt","r")
for line in f:
    line = line.strip().split("\t")
    sql = "INSERT INTO SimilarSequences VALUES ('%s','%s','%s','%s','%s','%s','%s','%s')" % (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7])
    cursor.execute(sql)
    db.commit()

