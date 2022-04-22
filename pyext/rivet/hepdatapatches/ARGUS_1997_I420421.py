import yoda,math
ARGUS_1997_I420421_hTemp3=None
ARGUS_1997_I420421_hTemp4=None
# add up the partial waves to get the total
def patch(path, ao):
    global ARGUS_1997_I420421_hTemp3,ARGUS_1997_I420421_hTemp4
    if "ARGUS_1997_I420421" in path:
        if "d04" in path :
            if ARGUS_1997_I420421_hTemp4 :
                for i in range(0,len(ao.points())) :
                    ARGUS_1997_I420421_hTemp4.points()[i].setY(ARGUS_1997_I420421_hTemp4.points()[i].y()+
                                                               ao.points()[i].y())
                    for (key,value) in ao.points()[i].errMap().items() :
                        value2=ARGUS_1997_I420421_hTemp4.points()[i].errMap()[key]
                        lower = math.sqrt(value[0]**2+value2[0]**2)
                        upper = math.sqrt(value[1]**2+value2[1]**2)
                        ARGUS_1997_I420421_hTemp4.points()[i].setYErrs((lower,upper),key)
                ARGUS_1997_I420421_hTemp4.updateTotalUncertainty()
                return[ao,ARGUS_1997_I420421_hTemp4]
            else :
                ARGUS_1997_I420421_hTemp4=yoda.Scatter2D()
                ARGUS_1997_I420421_hTemp4.setPath(ao.path().replace("y01","y03"))
                for p in ao.points():             
                    ARGUS_1997_I420421_hTemp4.addPoint(p.x(),p.y())
                    ARGUS_1997_I420421_hTemp4.points()[-1].setXErrs(p.xErrs())
                    for (key,value) in p.errMap().items() :
                        ARGUS_1997_I420421_hTemp4.points()[-1].setYErrs(value,key)
        elif "d03" in path :
            if ARGUS_1997_I420421_hTemp3 :
                ioff=0
                if "y04" in path : ioff=13
                for i in range(0,len(ao.points())) :
                    ARGUS_1997_I420421_hTemp3.points()[i+ioff].setY(ARGUS_1997_I420421_hTemp3.points()[i+ioff].y()+
                                                               ao.points()[i].y())
                    for (key,value) in ao.points()[i].errMap().items() :
                        value2=ARGUS_1997_I420421_hTemp3.points()[i+ioff].errMap()[key]
                        lower = math.sqrt(value[0]**2+value2[0]**2)
                        upper = math.sqrt(value[1]**2+value2[1]**2)
                        ARGUS_1997_I420421_hTemp3.points()[i+ioff].setYErrs((lower,upper),key)
                if "y04" in path :
                    ARGUS_1997_I420421_hTemp3.updateTotalUncertainty()
                    return[ao,ARGUS_1997_I420421_hTemp3]
            else :
                ARGUS_1997_I420421_hTemp3=yoda.Scatter2D()
                ARGUS_1997_I420421_hTemp3.setPath(ao.path().replace("y01","y05"))
                for p in ao.points():
                    ARGUS_1997_I420421_hTemp3.addPoint(p.x(),p.y())
                    ARGUS_1997_I420421_hTemp3.points()[-1].setXErrs(p.xErrs())
                    for (key,value) in p.errMap().items() :
                        ARGUS_1997_I420421_hTemp3.points()[-1].setYErrs(value,key)
    return ao
