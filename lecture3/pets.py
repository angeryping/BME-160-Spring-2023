class Pets:
    def __init__(self, name, color):
        self.color = color  # whatever color input in will be into the init
        self.name = name 
    
    def talk(self): #self refers to object it is working on right now
        print("My name is {} and i am {}".format(self.name, self.color))

def main():
    firstPet = Pets ('lucky', 'white')
    secondPet = Pets ('Charlie', 'brown')
    x = Pets ("Ralph",'purple')

    firstPet.talk()
    x.talk()
    secondPet.talk()
    
main()