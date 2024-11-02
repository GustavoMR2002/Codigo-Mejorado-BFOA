import blosum as bl

class evaluadorBlosum():
    
    def __init__(self):
        matrix = bl.BLOSUM(62)
        
        self.matrix = matrix
        
    def showMatrix(self):
        print(self.matrix)
        
    def getScore(self, A, B):
        score = self.matrix[A][B]
        #print("Score: ",score)
        return score
    
    
    pass




