'''
Il bit vector indicizza da 0. 
Presenta un 1 dove nel testo di trova un "|"
'''

class BV:
    def __init__(self, text):
        self.text = text
        
        self.length = len(text)
        
        self.bit_vector = []
        self.ones = []
        i = 1
        for c in self.text:
            if c == "|":
                self.bit_vector.append(1)
                self.ones.append(i)
            else:
                self.bit_vector.append(0)
            i += 1

    def rank(self, i):
        n = 0
        for elem in self.ones:
            if elem<=int(i):
                n+=1
            else:
                break
        return n
    
    def select(self, i):
        iterator = 0
        for elem in self.ones:
            iterator+=1
            if iterator == i:
                return elem

    def __str__(self):
        s = ""
        for elem in self.bit_vector:
            s += str(elem)
        return s

if __name__ == '__main__':
    pass
