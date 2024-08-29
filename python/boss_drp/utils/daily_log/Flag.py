
class Flag:
    def __init__(self,name,color,code,desc):
        self.name=name
        self.color=color
        self.code=code
        self.desc=desc
    def key(self):
        return(f"      <TR><TD><p style='color:{self.color};'>{self.name}</p><TD>{self.desc}</TR>\n")
        
incomplete = Flag('MAGENTA','magenta','#FF00FF','Reduction has yet to be started due to incomplete transfer')
stopped = Flag('RED','red','#FF0000','Stopped Reduction')
NoExp = Flag('MAROON','Maroon','#FFA500','Stopped Reduction for NO GOOD EXPOSURES')
Error_warn = Flag('ORANGE','DarkOrange','#FF8C00','Pipeline ran with errors/warnings')
running = Flag('YELLOW','Gold', '#FFD700','Pipeline is still running')
#NoRedux = Flag('BLUE','blue','#0000FF','No Reductions')
NoObs = Flag('BLUE','blue','#0000FF','No Observations')
NoRedux = Flag('TEAL','MediumTurquoise','#48D1CC','No Science Observations')
#NoObs = Flag('TEAL','Teal','#008080','No Observations')
NoIssues = Flag('GREEN','green','#008000','No issues')
Silent_warn = Flag('ORANGE','#FF8C00','#FF8C00','Pipeline ran with acceptable warnings')
