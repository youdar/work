def run():
    g = open('files_with_problems', "r").readlines()
    l1 = [x for x in g if 'NoneType' not in x]
    l2 = [x for x in g if 'NoneType' in x]
    #add_to_tested_files(l2)
    #clean_problem_files(l1)

    
def add_to_tested_files(data):
    data = [x[:4] + '::0::NO R-work information in pdb file!!!' for x in data]
    for x in data:
        print x    
    f = open('Collect_tested_files',"a")
    for d in data:
        f.write(d + '\n')    
    f.close()   
    
def clean_problem_files(data):
    #for x in data:
        #print x    

    g = open('files_with_problems',"w")
    for d in data:
        g.write(d)    
    g.close()    


if __name__=='__main__':
    run()