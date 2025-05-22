#*****************************************************************
# PROJECT 3: Von Der Malsburg - by Josh Goh 28 Apr 2007
#*****************************************************************

# Specify output file paths
outputfilepath=sprintf("%s/vdm_output_jg.txt",getwd())
outputfile<-file(outputfilepath,"a")

# Set project global variables
hexstartendL1=cbind(c(4,4,3,3,2,2,1,1,1,2,2,3,3,4,4),c(11,12,12,13,13,14,14,15,14,14,13,13,12,12,11))
hexstartendL2=cbind(c(2,2,1,2,2),c(4,5,5,5,4))
hexstartendstim=hexstartendL2
W=c(19,15); p=0.4; q=0.3; r=0.286; s=0.25; h=0.05; theta=1; nsect=9; order=c(1,1,2)


#-----------------------------------------------------------------
# READ FUNCTIONS INTO WORKSPACE
#-----------------------------------------------------------------

# Function 1: vdmulearning_2layer
#-----------------------------------------------------------------
# Trains a network according to the unsupervised learning network described in Von Der Malsburg (1973).
# Simulates visual cortex orientation cells. Hexagonal network arrangements of units is used to simulate equal
# contribution in surrounding units to each unit. This specific network is different from VDM in that it
# includes another layer that integrates orientation information from the original VDM network and that it
# also simulates several VDM networks representing retinotopic location in space.
# 
# Usage:
#
# vdmulearning_2layer(W,p,r,q,s,h,t,order,input,theta,hexstartendL1,hexstartendL2,nsect)
#
# W	        - initial weight matrix list with specified items:
#                 1 - Layer 1 Input-Excitatory weights (i-E)
#                 2 - Layer 1 Excitatory-Excitatory weights (E-E)
#                 3 - Layer 1 Excitatory-Inhibitory weights (E-I)
#                 4 - Layer 1 Inhibitory-Excitatory weights (I-E)
#                 5 - Layer 2 Input-Excitatory weights (i-E)
#                 6 - Layer 2 Excitatory-Excitatory weights (E-E)
#                 7 - Layer 2 Excitatory-Inhibitory weights (E-I)
#                 8 - Layer 2 Inhibitory-Excitatory weights (I-E)
#                 If no specified list, then set as vector with number of
#                 input units first, then the width of hexagonal layer 1 units,
#                 the the width of hexagonal layer 2 units.
# p             - strength of excitatory-excitatory weights
# r             - strength of excitatory-inhibitory weights
# q             - strength of inhibitory-excitatory weights
# s             - max initial strength of input-excitatory weight
# h             - growth rate
# t             - number of iterations
# order         - vector of orders of influence for cortical weights
#                 1 - E-E
#                 2 - E-I
#                 3 - I-E
# input         - list of input column vector(s) for each retinotopic location
# theta         - activation function threshold
# hexstartendL1 - vector of start (col1) and end (col2) units for each row of the
#                 considered layer 1 hexagon
# hexstartendL2 - vector of start (col1) and end (col2) units for each row of the
#                 considered layer 2 hexagon
# nsect         - number of layer 1 hexagons as input to layer 2

vdmulearning_2layer<-function(W,p,r,q,s,h,t,order,input,theta,hexstartendL1,hexstartendL2,nsect) {
		
	# Initialize weight connections
	
	WiE = c() # i-E weights
	WEE = c() # E-E weights
	WEI = c() # E-I weights
	WIE = c() # I-E weights
	ni  = c() # number of input units
	nE  = c() # number of excitatory units
	nI  = c() # number of inhibitory units
	
	if (!is.list(W)) {
		
		# Parameters
		# For general use
		size          = c()
		size[[1]]     = W[2]                                                                  # width of layer 1 units
		size[[2]]     = W[3]                                                                  # width of layer 2 units
		hexmat        = c()		
		hexmat[[1]]   = matrix(0,size[[1]][1],size[[1]][1])                                   # matrix of layer 1 hexagon indices
		hexmat[[2]]   = matrix(0,size[[2]][1],size[[2]][1])                                   # matrix of layer 2 hexagon indices
		hexmatno      = c()
		hexmatno[[1]] = t(matrix(seq(1,size[[1]][1]*size[[1]][1]),size[[1]][1],size[[1]][1])) # vector of layer 1 hexagon unit number
		hexmatno[[2]] = t(matrix(seq(1,size[[2]][1]*size[[2]][1]),size[[2]][1],size[[2]][1])) # vector of layer 2 hexagon unit number
		col1          = c()
		col1[[1]]     = hexstartendL1[,1]                                                     # start units for each layer 1 hexagonal row
		col2          = c()
		col2[[1]]     = hexstartendL1[,2]                                                     # end units for each layer 1 hexagonal row
		col1[[2]]     = hexstartendL2[,1]                                                     # start units for each layer 2 hexagonal row
		col2[[2]]     = hexstartendL2[,2]                                                     # end units for each layer 2 hexagonal row		
		ni[[1]]       = W[1]                                                                  # number of input units per layer 1 section
		unitsep       = 1                                                                     # set desired unit separation distance between hexagon units
		eesep         = order[1]                                                              # order of units of influence for E-E weights
		eisep         = order[2]                                                              # order of units of influence for E-I weights
		iesep         = order[3]                                                              # order of units of influence for I-E weights
		
		# For geometric hexagonal calculations
		hl       = 1            # length of hexagon edge
		h        = sin(pi/6)*hl # height of hexagon side hypotenuse
		w        = cos(pi/6)*hl # width of hexagon side hypotenuse

		# Convert hexmat indices to hexagonal corner indices in matrix space for both layers
		xycc=c()
		for (lay in 1:2) {
			xycctemp=c()
			for (hy in 1:dim(hexmat[[lay]])[1]) {
				for (hx in col1[[lay]][hy]:col2[[lay]][hy]) {
					if (hy%%2) {
						xcc=hx*2*w+w
						ycc=hy*(h+hl)
					} else {
						xcc=hx*2*w
						ycc=hy*(h+hl)
					}
					xycctemp=rbind(xycctemp,c(ycc,xcc))
				}
			}
			xycc[[lay]]=xycctemp
		} # end layer loop
		
        	# Define actual matrix distance between two hexagon centers
		hexrad=sqrt((xycc[[1]][1,1]-xycc[[1]][(unitsep+1),1])^2+(xycc[[1]][1,2]-xycc[[1]][(unitsep+1),2])^2)

		# Set excitatory, inhibitory connection weights for both layers
		for (lay in 1:2) {
			nE[[lay]]  = dim(xycc[[lay]])[1]         # Number of excitatory units
			nI[[lay]]  = nE[[lay]]                   # Number of inhibitory units
			WEE[[lay]] = c()
			WEI[[lay]] = c()
			WIE[[lay]] = c()
			WiE[[lay]] = c()
		} # end layer loop
		
		# Initialize layer 1 sections
		for (sect in 1:nsect) {
		
			WiE[[1]][[sect]] = matrix(runif(ni[[1]]*nE[[1]],0,s),ni[[1]],nE[[1]]) # i-E
			WEE[[1]][[sect]] = matrix(0,nE[[1]],nE[[1]])                          # E-E weights
		        WEI[[1]][[sect]] = WEE[[1]][[sect]]                                   # E-I weights
		        WIE[[1]][[sect]] = WEE[[1]][[sect]]                                   # I-E weights
		
			for (u1 in 1:nE[[1]]) {
				for (u2 in 1:nE[[1]]) {
					dy=xycc[[1]][u2,1]-xycc[[1]][u1,1]
					dx=xycc[[1]][u2,2]-xycc[[1]][u1,2]

					# E-E
					if (round(sqrt(dx^2+dy^2))==round(eesep*hexrad)) {
						WEE[[1]][[sect]][u1,u2]=p
					}

					# E-I
					if (round(sqrt(dx^2+dy^2))<=round(eisep*hexrad)) {
						WEI[[1]][[sect]][u1,u2]=r
					}

					# I-E
					if (round(sqrt(dx^2+dy^2))==round(iesep*hexrad)) {
						WIE[[1]][[sect]][u1,u2]=q
					}

				}
			}
			dy=c();dx=c()
			
		} # End section loop
		
		# Initialize layer 2
		ni[[2]]  = nsect*nE[[1]]
		WiE[[2]] = matrix(runif(ni[[2]]*nE[[2]],0,s),ni[[2]],nE[[2]]) # i-E
		WEE[[2]] = matrix(0,nE[[2]],nE[[2]])                          # E-E weights
	        WEI[[2]] = WEE[[2]]                                           # E-I weights
	        WIE[[2]] = WEE[[2]]                                           # I-E weights
		
		for (u1 in 1:nE[[2]]) {
				for (u2 in 1:nE[[2]]) {
					dy=xycc[[2]][u2,1]-xycc[[2]][u1,1]
					dx=xycc[[2]][u2,2]-xycc[[2]][u1,2]

					# E-E
					if (round(sqrt(dx^2+dy^2))==round(eesep*hexrad)) {
						WEE[[2]][u1,u2]=p
					}

					# E-I
					if (round(sqrt(dx^2+dy^2))<=round(eisep*hexrad)) {
						WEI[[2]][u1,u2]=r
					}

					# I-E
					if (round(sqrt(dx^2+dy^2))==round(iesep*hexrad)) {
						WIE[[2]][u1,u2]=q
					}

				}
			}
		
	} else {
		
		WiE[[1]] = W[[1]]
		WEE[[1]] = W[[2]]
		WEI[[1]] = W[[3]]
		WIE[[1]] = W[[4]]
		WiE[[2]] = W[[5]]
		WEE[[2]] = W[[6]]
		WEI[[2]] = W[[7]]
		WIE[[2]] = W[[8]]
		ni[[1]]  = dim(WiE[[1]])[1]
		nE[[1]]  = dim(WiE[[1]])[2]
		nI[[1]]  = nE[[1]]
		ni[[2]]  = dim(WiE[[2]])[1]
		nE[[2]]  = dim(WiE[[2]])[2]
		nI[[2]]  = nE[[2]]
		
	} # End weight initialization
	
	
	# Train network	
	# Initialize Parameters
	asEt      = c()            # excitatory states over runs
	asIt	  = c()            # inhibitory states over runs
	asEdt     = c()            # excitatory states over every iteration
	asIdt     = c()            # inhibitory states over every iteration
	sigE      = c()            # excitatory output signals
	sigI      = c()            # inhibitory output signals	
	WiEt      = c()            # list of input-excitatory unit weights
	WEEt      = c()            # list of excitatory-excitatory unit weights
	WEIt      = c()            # list of excitatory-inhibitory unit weights
	WIEt      = c()            # list of inhibitory-excitatory unit weights
	maxast    = c()            # max unit difference between two iterations
	ip        = c()
	excite    = c()
	inhib     = c()
	sumsk     = c()
	WiEp      = c()
	asE       = c()
	asI       = c()
	WiEtlay=c()
	WEEtlay=c()
	WEItlay=c()
	WIEtlay=c()
	asEtlay=c()
	asItlay=c()
	WiEtsect=c()
	WEEtsect=c()
	WEItsect=c()
	WIEtsect=c()
	asEtsect=c()
	asItsect=c()
	
	for (sect in 1:nsect) {
		sigE[[1]][[sect]]   = matrix(0,nE[[1]],1) # excitatory output signals layer 1
		sigI[[1]][[sect]]   = matrix(0,nI[[1]],1) # inhibitory output signals layer 1
		asEdt[[1]][[sect]]  = matrix(0,nE[[1]],1) # excitatory states over every iteration layer 2
		asIdt[[1]][[sect]]  = matrix(0,nI[[1]],1) # inhibitory states over every iteration layer 2
		maxast[[1]][[sect]] = matrix(0,nE[[1]],1) # max unit difference between two iterations
		ip[[1]][[sect]]     = matrix(0,nE[[1]],1)
		excite[[1]][[sect]] = matrix(0,nE[[1]],1)
		inhib[[1]][[sect]]  = matrix(0,nE[[1]],1)
	}
	
	sigE[[2]]   = matrix(0,nE[[2]],1) # excitatory output signals layer 2
	sigI[[2]]   = matrix(0,nI[[2]],1) # inhibitory output signals layer 2
	asEdt[[2]]  = matrix(0,nE[[2]],1) # excitatory states over every iteration layer 2
	asIdt[[2]]  = matrix(0,nI[[2]],1) # inhibitory states over every iteration layer 2
	maxast[[2]] = matrix(0,nE[[2]],1) # max unit difference between two iterations
	ip[[2]]     = matrix(0,nE[[2]],1)
	excite[[2]] = matrix(0,nE[[2]],1)
	inhib[[2]]  = matrix(0,nE[[2]],1)
	
	# Start runs	
	for (run in 1:t) {
	
		asEp = c() # excitatory states for all patterns in a run
		asIp = c() # inhibitory states for all patterns in a run		
		
		# Set layer start states
		for (lay in 1:2) {
			if (lay==1) {
				for (sect in 1:nsect) {
					# Normalize weights
					sumsk[[lay]][[sect]] = colSums(WiE[[lay]][[sect]])
					WiE[[lay]][[sect]]   = WiE[[lay]][[sect]]*matrix(ni[[lay]]*s/2/sumsk[[lay]][[sect]],ni[[lay]],nE[[lay]],byrow=TRUE)
					WiEp[[lay]][[sect]]  = WiE[[lay]][[sect]]*0
					
					asEp[[lay]][[sect]] = matrix(0,nE[[1]],1) # excitatory states for all patterns in a run
					asIp[[lay]][[sect]] = matrix(0,nI[[1]],1) # inhibitory states for all patterns in a run
					
				} # end section loop
			} else {
				# Normalize weights
				sumsk[[lay]] = colSums(WiE[[lay]])
				WiE[[lay]]   = WiE[[lay]]*matrix(ni[[lay]]*s/2/sumsk[[lay]],ni[[lay]],nE[[lay]],byrow=TRUE)
				WiEp[[lay]]  = WiE[[lay]]*0
				
				asEp[[lay]] = matrix(0,nE[[2]],1) # excitatory states for all patterns in a run
				asIp[[lay]] = matrix(0,nI[[2]],1) # inhibitory states for all patterns in a run				
			}						
		} # end layer loop
		
		# Loop through input patterns
		for (pat in 1:dim(input[[1]])[2]) {
			
			for (sect in 1:nsect) {
				asE[[1]][[sect]] = matrix(0,nE[[1]],1) # layer 1 units excitatory states
				asI[[1]][[sect]] = matrix(0,nI[[1]],1) # layer 1 units inhibitory states
			} # end section loop
			
			asE[[2]] = matrix(0,nE[[2]],1) # layer 2 units excitatory states
			asI[[2]] = matrix(0,nI[[2]],1) # layer 2 units inhibitory states
												
			# Begin stabilizing network based on input
			for (count in 1:20) {
			
				L2input=c()
			
				# Layer 1
				for (sect in 1:nsect) {
				
					inputp = input[[sect]][,pat]

					# Log activation states
					asEdt[[1]][[sect]] = cbind(asEdt[[1]][[sect]],asE[[1]][[sect]])
					asIdt[[1]][[sect]] = cbind(asIdt[[1]][[sect]],asI[[1]][[sect]])

					# Compute excitatory output signals
					for (j in 1:nE[[1]]) {		
						if (asE[[1]][[sect]][j]>theta) {
							sigE[[1]][[sect]][j]=asE[[1]][[sect]][j]-theta
						}
						else {
							sigE[[1]][[sect]][j]=0
						}		
					}

					# Compute inhibitory output signals
					asI[[1]][[sect]]=t(WEI[[1]][[sect]])%*%sigE[[1]][[sect]]
					for (j in 1:nI[[1]]) {
						if (asI[[1]][[sect]][j]>theta) {
							sigI[[1]][[sect]][j]=asI[[1]][[sect]][j]-theta
						}
						else {
							sigI[[1]][[sect]][j]=0
						}
					}

					# Compute excitatory activation states
					ip[[1]][[sect]]     = cbind(ip[[1]][[sect]],(t(WiE[[1]][[sect]])%*%inputp))
					excite[[1]][[sect]] = cbind(excite[[1]][[sect]],(t(WEE[[1]][[sect]])%*%sigE[[1]][[sect]]))
					inhib[[1]][[sect]]  = cbind(inhib[[1]][[sect]],(t(WIE[[1]][[sect]])%*%sigI[[1]][[sect]]))
					asEnew              = (t(WiE[[1]][[sect]])%*%inputp)+(t(WEE[[1]][[sect]])%*%sigE[[1]][[sect]])-(t(WIE[[1]][[sect]])%*%sigI[[1]][[sect]])					
					maxast[[1]][[sect]] = cbind(maxast[[1]][[sect]],(asEnew-asE[[1]][[sect]]))
					asE[[1]][[sect]]    = asEnew
					
					# Compute input to layer 2
					L2input=rbind(L2input,asE[[1]][[sect]])

				} # End section loop
				
				
			
				# Layer 2
				# Binarize input from layer 1
				for (i in 1:length(L2input)) {
					if (L2input[i]>theta) {
						L2input[i]=1
					} else {
						L2input[i]=0
					}
				}
				
				
				# Log activation states
				asEdt[[2]]  = cbind(asEdt[[2]],asE[[2]])
				asIdt[[2]]  = cbind(asIdt[[2]],asI[[2]])

				# Compute excitatory output signals
				for (j in 1:nE[[2]]) {		
					if (asE[[2]][j]>theta) {
						sigE[[2]][j]=asE[[2]][j]-theta
					}
					else {
						sigE[[2]][j]=0
					}		
				}

				# Compute inhibitory output signals
				asI[[2]]=t(WEI[[2]])%*%sigE[[2]]
				for (j in 1:nI[[2]]) {
					if (asI[[2]][j]>theta) {
						sigI[[2]][j]=asI[[2]][j]-theta
					}
					else {
						sigI[[2]][j]=0
					}
				}

				# Compute excitatory activation states
				ip[[2]]      = cbind(ip[[2]],(t(WiE[[2]])%*%L2input))
				excite[[2]]  = cbind(excite[[2]],(t(WEE[[2]])%*%sigE[[2]]))
				inhib[[2]]   = cbind(inhib[[2]],(t(WIE[[2]])%*%sigI[[2]]))
				asEnew      = (t(WiE[[2]])%*%L2input)+(t(WEE[[2]])%*%sigE[[2]])-(t(WIE[[2]])%*%sigI[[2]])
				maxast[[2]]  = cbind(maxast[[2]],(asEnew-asE[[2]]))
				asE[[2]]    = asEnew
			
			} # End stabilizing loop
			
			# Updates
			for (lay in 1:2) {
				if (lay==1) {
					for (sect in 1:nsect) {
						# Log activation states
						asEdt[[lay]][[sect]]  = cbind(asEdt[[lay]][[sect]],asE[[lay]][[sect]])
						asIdt[[lay]][[sect]]  = cbind(asIdt[[lay]][[sect]],asI[[lay]][[sect]])
						
						# Update weights and unit activation states for each pattern
						asEp[[lay]][[sect]] = cbind(asEp[[lay]][[sect]],asE[[lay]][[sect]])
						asIp[[lay]][[sect]] = cbind(asIp[[lay]][[sect]],asI[[lay]][[sect]])			
						WiEp[[lay]][[sect]] = WiEp[[lay]][[sect]]+h*inputp%*%t(sigE[[lay]][[sect]])
						
					} # end section loop
				} else {
					# Log activation states
					asEdt[[lay]]  = cbind(asEdt[[lay]],asE[[lay]])
					asIdt[[lay]]  = cbind(asIdt[[lay]],asI[[lay]])

					# Update weights and unit activation states for each pattern
					asEp[[lay]] = cbind(asEp[[lay]],asE[[lay]])
					asIp[[lay]] = cbind(asIp[[lay]],asI[[lay]])			
					WiEp[[lay]] = WiEp[[lay]]+h*L2input%*%t(sigE[[lay]])
					
				}						
			} # end layer loop
						
		} # End input pattern loop
		
		# Updates
		for (lay in 1:2) {
			if (lay==1) {
				for (sect in 1:nsect) {
					WiE[[lay]][[sect]] = WiE[[lay]][[sect]]+WiEp[[lay]][[sect]]
					WiEtsect[[sect]] = WiE[[lay]][[sect]]
					WEEtsect[[sect]] = WEE[[lay]][[sect]]
					WEItsect[[sect]] = WEI[[lay]][[sect]]
					WIEtsect[[sect]] = WIE[[lay]][[sect]]
					asEtsect[[sect]] = asEp[[lay]][[sect]]
					asItsect[[sect]] = asIp[[lay]][[sect]]
				} # end section loop
				
				WiEtlay[[lay]]=WiEtsect
				WEEtlay[[lay]]=WEEtsect
				WEItlay[[lay]]=WEItsect
				WIEtlay[[lay]]=WIEtsect
				asEtlay[[lay]]=asEtsect
				asItlay[[lay]]=asItsect
				
			} else {
				WiE[[lay]] = WiE[[lay]]+WiEp[[lay]]
				WiEtlay[[lay]] = WiE[[lay]]
				WEEtlay[[lay]] = WEE[[lay]]
				WEItlay[[lay]] = WEI[[lay]]
				WIEtlay[[lay]] = WIE[[lay]]
				asEtlay[[lay]] = asEp[[lay]]
				asItlay[[lay]] = asIp[[lay]]
			}
		} # end layer loop
		
		WiEt[[run]]=WiEtlay
		WEEt[[run]]=WEEtlay
		WEIt[[run]]=WEItlay
		WIEt[[run]]=WIEtlay
		asEt[[run]]=asEtlay
		asIt[[run]]=asItlay
								
	} # End run loop
	
	# Remove extra columns in variables
	for (time in 1:run) {
		for (lay in 1:2) {
			if (lay==1) {
				for (sect in 1:nsect) {
					asEt[[time]][[lay]][[sect]]=asEt[[time]][[lay]][[sect]][,-1]
					#asIt[[time]][[lay]][[sect]]=asIt[[time]][[lay]][[sect]][,-1]
				}
			} else {
				asEt[[time]][[lay]]=asEt[[time]][[lay]][,-1]
			}
		}
	}
	
	
	# Update global variables
	assign("WiEt",WiEt,envir=.GlobalEnv)
	assign("WEEt",WEEt,envir=.GlobalEnv)
	assign("WEIt",WEIt,envir=.GlobalEnv)
	assign("WIEt",WIEt,envir=.GlobalEnv)
	assign("asEt",asEt,envir=.GlobalEnv)
	#assign("asIt",asIt,envir=.GlobalEnv)
	#assign("maxast",maxast,envir=.GlobalEnv)
	#assign("asEdt",asEdt,envir=.GlobalEnv)
	#assign("asIdt",asIdt,envir=.GlobalEnv)
	#assign("ip",ip,envir=.GlobalEnv)
	#assign("excite",excite,envir=.GlobalEnv)
	#assign("inhib",inhib,envir=.GlobalEnv)
	
		
} # End function 1


# Function 2: convert VDM activation to hexagon plot
#-----------------------------------------------------------------

vdmhexplot<-function(hexstartend,asE,theta,hexsize) {

	# Parameters
	hl       = hexsize             # length of hexagon edge
	h        = sin(pi/6)*hl        # height of hexagon side hypotenuse
	w        = cos(pi/6)*hl        # width of hexagon side hypotenuse
	hexw     = dim(hexstartend)[1] # hexagonal width of units
	xycc     = c()
	x1       = c()
	x2       = c()
	y1       = c()
	y2       = c()
	u        = 1
	i        = 1
	j        = 1

	# Compute hexagon state
	for (hy in hexw:1) {
		for (hx in hexstartend[hy,1]:hexstartend[hy,2]) {
			if (asE[u]>=theta) {
				if (hy%%2) {
					x1[i]=hx*2*w+w		
				} else {
					x1[i]=hx*2*w
				}
				y1[i]=hy*(h+hl)
				i=i+1
			} else {
				if (hy%%2) {
					x2[j]=hx*2*w+w
				} else {
					x2[j]=hx*2*w
				}
				y2[j]=hy*(h+hl)
				j=j+1
			}
			u=u+1
		}
	}
	
	# Assign variable to global environment
	assign("asEhexcoordsON",cbind(y1,x1),envir=.GlobalEnv)
	assign("asEhexcoordsOFF",cbind(y2,x2),envir=.GlobalEnv)
	
	# Plotting
	cc=c(x1,x2,y1,y2)
	plot(min(cc):max(cc),min(cc):max(cc),xlab="",ylab="",type="n",axes=FALSE)
	points(x2,y2,pch=19,col="black",cex=0.5)
	if (!is.null(asEhexcoordsON)) {
	points(x1,y1,pch=19,col="black",cex=1.5)
	}
	
} # end function 2

# Function 3: Train Unsupervised Learning Network basic
#
# Trains a network according to the unsupervised learning network described in Von Der Malsburg (1973). 
# Simulates visual cortex orientation cells. Hexagonal network arrangements of units is used so simulate
# equal contribution in surrounding units to each unit.
# 
# Usage:
#
# vdmulearning(W,p,r,q,s,h,t,order,input,theta,hexstartend)
#
# W	      - initial weight matrix list with specified items:
#               1 - Input-Excitatory weights (i-E)
#               2 - Excitatory-Excitatory weights (E-E)
#               3 - Excitatory-Inhibitory weights (E-I)
#               4 - Inhibitory-Excitatory weights (I-E)
#               If no specified list, then set as vector with number of
#               input units first, then the width of hexagonal layer 2 units.
# p           - strength of excitatory-excitatory weights
# r           - strength of excitatory-inhibitory weights
# q           - strength of inhibitory-excitatory weights
# s           - max initial strength of input-excitatory weight
# h           - growth rate
# t           - number of iterations
# order       - vector of orders of influence for cortical weights
#               1 - E-E
#               2 - E-I
#               3 - I-E
# input       - input column vector(s)
# theta       - activation function threshold
# hexstartend - vector of start (col1) and end (col2) units for each row of the considered hexagon

vdmulearning<-function(W,p,r,q,s,h,t,order,input,theta,hexstartend) {
		
	# Initialize weight connections
	if (!is.list(W)) {
		
		# Parameters
		# For general use
		size     = W[2]                                  # no. of layer 2 units
		eesep    = order[1]                              # order of units of influence for E-E weights
		eisep    = order[2]                              # order of units of influence for E-I weights
		iesep    = order[3]                              # order of units of influence for I-E weights
		hexmat   = matrix(0,size,size)                   # matrix of layer 2 hexagon indices
		hexmatno = t(matrix(seq(1,size*size),size,size)) # vector of layer 2 hexagon unit number
		unitsep  = 1                                     # set desired unit separation distance between hexagon units
		col1     = hexstartend[,1]                       # start units for each hexagonal row
		col2     = hexstartend[,2]                       # end units for each hexagonal row
		ni       = W[1]                                  # number of input units
		
		
		# For geometric hexagonal calculations
		hl       = 1                                     # length of hexagon edge
		h        = sin(pi/6)*hl                          # height of hexagon side hypotenuse
		w        = cos(pi/6)*hl                          # width of hexagon side hypotenuse

		# Convert hexmat indices to hexagonal corner indices in matrix space
		xycc=c()
		for (hy in 1:dim(hexmat)[1]) {
			for (hx in col1[hy]:col2[hy]) {
				if (hy%%2) {
					xcc=hx*2*w+w
					ycc=hy*(h+hl)
				} else {
					xcc=hx*2*w
					ycc=hy*(h+hl)
				}
				xycc=rbind(xycc,c(ycc,xcc))
			}
		}
		
                # Define actual matrix distance between two hexagon centers
		hexrad=sqrt((xycc[1,1]-xycc[(unitsep+1),1])^2+(xycc[1,2]-xycc[(unitsep+1),2])^2)

		# Set excitatory, inhibitory connection weights
		nE  = dim(xycc)[1]         # Number of excitatory units
		nI  = nE                   # Number of inhibitory units
		WEE = matrix(0,nE,nE)      # E-E weights
		WEI = WEE                  # E-I weights
		WIE = WEE                  # I-E weights
		
		# i-E
		WiE=matrix(runif(ni*nE,0,s),ni,nE)
		
		for (u1 in 1:nE) {
			for (u2 in 1:nE) {
				dy=xycc[u2,1]-xycc[u1,1]
				dx=xycc[u2,2]-xycc[u1,2]

				# E-E
				if (round(sqrt(dx^2+dy^2))==round(eesep*hexrad)) {
					WEE[u1,u2]=p
				}

				# E-I
				if (round(sqrt(dx^2+dy^2))<=round(eisep*hexrad)) {
					WEI[u1,u2]=r
				}

				# I-E
				if (round(sqrt(dx^2+dy^2))==round(iesep*hexrad)) {
					WIE[u1,u2]=q
				}

			}
		}
		
	} else {
		WiE=W[[1]]
		WEE=W[[2]]
		WEI=W[[3]]
		WIE=W[[4]]
		
		# Set parameters
		ni=dim(WiE)[1]	# Number of input units
		nE=dim(WiE)[2]	# Number of excitatory units
		nI=nE		# Number of inhibitory units
	}
	
	
	# Train network	
	# Initialize Parameters
	asEt      = c()            # excitatory states over runs
	asIt	  = c()            # inhibitory states over runs
	asEdt     = c()            # excitatory states over every iteration
	asIdt     = c()            # inhibitory states over every iteration
	sigE      = c()            # excitatory output signals
	sigI      = c()            # inhibitory output signals
	WiEt      = c()            # list of input-excitatory unit weights
	WEEt      = c()            # list of excitatory-excitatory unit weights
	WEIt      = c()            # list of excitatory-inhibitory unit weights
	WIEt      = c()            # list of inhibitory-excitatory unit weights
	maxast    = c()            # max unit difference between two iterations
	ip        = c()
	excite    = c()
	inhib     = c()
	
	# Start runs	
	for (run in 1:t) {
	
		asEp = c() # excitatory states for all patterns in a run
		asIp = c() # inhibitory states for all patterns in a run
		
		# Normalize weights
		sumsk = colSums(WiE)
		WiE   = WiE*matrix(ni*s/2/sumsk,ni,nE,byrow=TRUE)
		WiEp  = WiE*0
		
		# Loop through input patterns
		for (pat in 1:dim(input)[2]) {
			
			inputp = input[,pat]
			asE    = matrix(0,nE,1) # excitatory states
			asI    = matrix(0,nI,1) # inhibitory states			
			
			# Begin stabilizing network based on input		
			for (count in 1:20) {
				
				# Log activation states
				asEdt  = cbind(asEdt,asE)
				asIdt  = cbind(asIdt,asI)

				# Compute excitatory output signals
				for (j in 1:nE) {		
					if (asE[j]>theta) {
						sigE[j]=asE[j]-theta
					}
					else {
						sigE[j]=0
					}		
				}

				# Compute inhibitory output signals
				asI=t(WEI)%*%sigE
				for (j in 1:nI) {
					if (asI[j]>theta) {
						sigI[j]=asI[j]-theta
					}
					else {
						sigI[j]=0
					}
				}

				# Compute excitatory activation states
				ip     = cbind(ip,(t(WiE)%*%inputp))
				excite = cbind(excite,(t(WEE)%*%sigE))
				inhib  = cbind(inhib,(t(WIE)%*%sigI))
				asEnew = (t(WiE)%*%inputp)+(t(WEE)%*%sigE)-(t(WIE)%*%sigI) 
				maxast = cbind(maxast,(asEnew-asE))
				asE    = asEnew				
			
			} # End stabilizing loop
			
			# Log activation states
			asEdt  = cbind(asEdt,asE)
			asIdt  = cbind(asIdt,asI)
			
			# Update weights and unit activation states for each pattern
			asEp = cbind(asEp,asE)
			asIp = cbind(asIp,asI)			
			WiEp = WiEp+h*inputp%*%t(sigE)
			
		} # End input pattern loop
		
		# Update weights
		WiE = WiE+WiEp
		
		# Update variables over runs
		WiEt[[run]] = WiE
		WEEt[[run]] = WEE
		WEIt[[run]] = WEI
		WIEt[[run]] = WIE
		asEt[[run]] = asEp
		asIt[[run]] = asIp
		assign("WiEt",WiEt,envir=.GlobalEnv)
		assign("WEEt",WEEt,envir=.GlobalEnv)
		assign("WEIt",WEIt,envir=.GlobalEnv)
		assign("WIEt",WIEt,envir=.GlobalEnv)
		assign("asEt",asEt,envir=.GlobalEnv)
		assign("asIt",asIt,envir=.GlobalEnv)
		assign("maxast",maxast,envir=.GlobalEnv)
		assign("asEdt",asEdt,envir=.GlobalEnv)
		assign("asIdt",asIdt,envir=.GlobalEnv)
		assign("ip",ip,envir=.GlobalEnv)
		assign("excite",excite,envir=.GlobalEnv)
		assign("inhib",inhib,envir=.GlobalEnv)
		
			
	} # End run loop
		
} # End function 3


#-----------------------------------------------------------------
# BEGIN WORKSPACE COMPUTATIONS FOR PROJECT SIMULATIONS
#-----------------------------------------------------------------

# CREATE SIMULATION INPUT PATTERNS
#-----------------------------------------------------------------
# BASE ORIENTATION INPUT
input=            matrix(c(0,1,0,0,1,1,0,0,0,1,0,0,0,1,1,0,0,1,0),19,1)
input=cbind(input,matrix(c(0,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0),19,1))
input=cbind(input,matrix(c(0,0,1,0,0,1,1,0,0,1,0,0,1,1,0,0,1,0,0),19,1))
input=cbind(input,matrix(c(0,0,0,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,0),19,1))
input=cbind(input,matrix(c(0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0),19,1))
input=cbind(input,matrix(c(0,0,0,1,0,0,0,1,1,1,1,1,0,0,0,1,0,0,0),19,1))
input=cbind(input,matrix(c(0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,0,0,0),19,1))
input=cbind(input,matrix(c(1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,0,0,1),19,1))
input=cbind(input,matrix(c(1,1,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,1),19,1))

# Include retintopy into input
rinput=c()
for (l in 1:9) {
	rinput[[l]]=input
}

# Vertical only
vertinput=cbind(input[,2],input[,9])

# Include retintopy into input
rvinput=c()
for (l in 1:9) {
	rvinput[[l]]=vertinput
}

# Horizontal only
horzinput=cbind(input[,4],input[,5])

# Include retintopy into input
rhinput=c()
for (l in 1:9) {
	rhinput[[l]]=horzinput
}

# Print base input figure
png(file="base_input_fig.png")
layout(matrix(c(1:9),3,3,byrow=TRUE));layout.show(9)
for (i in 1:9) {vdmhexplot(hexstartendstim,input[,i],1,1)}
dev.off()

# SQUARE INPUTS (retinotopic)
squareinput=c()

# Square 1
squareinput[[1]]=input[,5]
squareinput[[2]]=input[,1]
squareinput[[3]]=input[,1]
squareinput[[4]]=matrix(0,19,1)
squareinput[[5]]=input[,6]
for (i in 6:9) {
	squareinput[[i]]=matrix(0,19,1)
}

# Square 2
squareinput[[1]]=cbind(squareinput[[1]],matrix(0,19,1))
squareinput[[2]]=cbind(squareinput[[2]],input[,5])
squareinput[[3]]=cbind(squareinput[[3]],matrix(0,19,1))
squareinput[[4]]=cbind(squareinput[[4]],input[,1])
squareinput[[5]]=cbind(squareinput[[5]],input[,1])
squareinput[[6]]=cbind(squareinput[[6]],matrix(0,19,1))
squareinput[[7]]=cbind(squareinput[[7]],input[,6])
for (i in 8:9) {
	squareinput[[i]]=cbind(squareinput[[i]],matrix(0,19,1))
}

# Square 3
squareinput[[1]]=cbind(squareinput[[1]],matrix(0,19,1))
squareinput[[2]]=cbind(squareinput[[2]],matrix(0,19,1))
squareinput[[3]]=cbind(squareinput[[3]],input[,5])
squareinput[[4]]=cbind(squareinput[[4]],matrix(0,19,1))
squareinput[[5]]=cbind(squareinput[[5]],input[,1])
squareinput[[6]]=cbind(squareinput[[6]],input[,1])
squareinput[[7]]=cbind(squareinput[[7]],matrix(0,19,1))
squareinput[[8]]=cbind(squareinput[[8]],input[,6])
squareinput[[9]]=cbind(squareinput[[9]],matrix(0,19,1))

# Square 4
for (i in 1:4) {
	squareinput[[i]]=cbind(squareinput[[i]],matrix(0,19,1))
}
squareinput[[5]]=cbind(squareinput[[5]],input[,5])
squareinput[[6]]=cbind(squareinput[[6]],matrix(0,19,1))
squareinput[[7]]=cbind(squareinput[[7]],input[,1])
squareinput[[8]]=cbind(squareinput[[8]],input[,1])
squareinput[[9]]=cbind(squareinput[[9]],input[,6])

# Print square input figure
for (l in 1:4) {

png(file=sprintf("square_%.0f_input_fig.png",l),width=800,height=800)
	layout(matrix(c(0,0,1,0,0,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,0,0,9,0,0),5,5,byrow=TRUE));layout.show(9)
	for (i in 1:9) {vdmhexplot(hexstartendstim,squareinput[[i]][,l],1,1)}
	dev.off()
} # end square location loop

# CROSS INPUTS
crossinput=c()

# Cross 1
crossinput[[1]]=input[,1]
crossinput[[2]]=input[,5]
crossinput[[3]]=input[,6]
crossinput[[4]]=matrix(0,19,1)
crossinput[[5]]=input[,1]
for (i in 6:9) {
	crossinput[[i]]=matrix(0,19,1)
}

# Cross 2
crossinput[[1]]=cbind(crossinput[[1]],matrix(0,19,1))
crossinput[[2]]=cbind(crossinput[[2]],input[,1])
crossinput[[3]]=cbind(crossinput[[3]],matrix(0,19,1))
crossinput[[4]]=cbind(crossinput[[4]],input[,5])
crossinput[[5]]=cbind(crossinput[[5]],input[,6])
crossinput[[6]]=cbind(crossinput[[6]],matrix(0,19,1))
crossinput[[7]]=cbind(crossinput[[7]],input[,1])
for (i in 8:9) {
	crossinput[[i]]=cbind(crossinput[[i]],matrix(0,19,1))
}

# Cross 3
crossinput[[1]]=cbind(crossinput[[1]],matrix(0,19,1))
crossinput[[2]]=cbind(crossinput[[2]],matrix(0,19,1))
crossinput[[3]]=cbind(crossinput[[3]],input[,1])
crossinput[[4]]=cbind(crossinput[[4]],matrix(0,19,1))
crossinput[[5]]=cbind(crossinput[[5]],input[,5])
crossinput[[6]]=cbind(crossinput[[6]],input[,6])
crossinput[[7]]=cbind(crossinput[[7]],matrix(0,19,1))
crossinput[[8]]=cbind(crossinput[[8]],input[,1])
crossinput[[9]]=cbind(crossinput[[9]],matrix(0,19,1))

# Cross 4
for (i in 1:4) {
	crossinput[[i]]=cbind(crossinput[[i]],matrix(0,19,1))
}
crossinput[[5]]=cbind(crossinput[[5]],input[,1])
crossinput[[6]]=cbind(crossinput[[6]],matrix(0,19,1))
crossinput[[7]]=cbind(crossinput[[7]],input[,5])
crossinput[[8]]=cbind(crossinput[[8]],input[,6])
crossinput[[9]]=cbind(crossinput[[9]],input[,1])

# Print cross input figure
for (l in 1:4) {
	png(file=sprintf("cross_%.0f_input_fig.png",l),width=800,height=800)
	layout(matrix(c(0,0,1,0,0,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,0,0,9,0,0),5,5,byrow=TRUE));layout.show(9)
	for (i in 1:9) {vdmhexplot(hexstartendstim,crossinput[[i]][,l],1,1)}
	dev.off()
} # end cross retinotopic location loop


# Combine square and cross inputs
scinput=c()
for (pos in 1:9) {
	scinput[[pos]]=cbind(squareinput[[pos]],crossinput[[pos]])
}


# PART 1: REPLICATE VDM ORIENTATION COLUMNS
#-----------------------------------------------------------------
vdmulearning(W,p,r,q,s,h,100,order,input,theta,hexstartendL1)

# Print layer 1 state figure at 1st and 100th run
for (k in c(1,100)) {
	for (l in 1:9) {
		png(file=sprintf("input_%.0f_layer1_r%.0f_fig.png",l,k))
		vdmhexplot(hexstartendL1,asEt[[k]][,l],1,1)
		dev.off()
	} 
}

# Compute unit modal classes
zeromod=matrix(0,1,3)
onemod=zeromod
polymod=zeromod
kc=1
for (k in c(1,20,100)) {
	for (u in 1:169) {
		actstate=c(1,asEt[[k]][u,],1)
		modcount=0
		for (pat in 2:10) {
			if (((actstate[pat]-actstate[pat-1])>0.5) & ((actstate[pat]-actstate[pat+1])>0.5) & (actstate[pat]>1)) {
				modcount=modcount+1
			}
		}
		if (modcount==0) {
			zeromod[1,kc]=zeromod[1,kc]+1
		}
		if (modcount==1) {
			onemod[1,kc]=onemod[1,kc]+1
		}
		if (modcount>1) {
			polymod[1,kc]=polymod[1,kc]+1
		}
	}
	kc=kc+1
}


# Print to output
cat("PART 1: REPLICATING VDM ORIENTATION COLUMNS IN LAYER 1\n",file=outputfilepath,append=TRUE)
cat("No. of units with 0, 1, >1 modes (rows) at 1, 20, and 100 time steps (cols)\n",file=outputfilepath,append=TRUE)
write.table(zeromod,sep="\t",file=outputfilepath,append=TRUE,col.names=FALSE)
cat("\n",file=outputfilepath,append=TRUE)
write.table(onemod,sep="\t",file=outputfilepath,append=TRUE,col.names=FALSE)
cat("\n",file=outputfilepath,append=TRUE)
write.table(polymod,sep="\t",file=outputfilepath,append=TRUE,col.names=FALSE)
cat("\n",file=outputfilepath,append=TRUE)	


# PART 2a: ANISOTROPY
#-----------------------------------------------------------------
vdmulearning(W,p,r,q,s,h,20,order,vertinput,theta,hexstartendL1)

# Print layer 1 state figure at 1st and 50th run

ent=matrix(0,2,2) # initiate entropy variable
kc=1

for (k in c(1,20)) {
	for (l in 1:2) {
		png(file=sprintf("vertinput_%.0f_layer1_r%.0f_fig.png",l,k))
		vdmhexplot(hexstartendL1,asEt[[k]][,l],1,1)
			
		# Compute excitatory activation for entropy
		ent[kc,l]=ent[kc,l]+sum(asEt[[k]])

		dev.off()
	}
	kc=kc+1
}

# Calculate mean entropy
vertment=ent/169

# Print to output
cat("PART 2: ANISOTROPY\n",file=outputfilepath,append=TRUE)
cat("Mean entropy at 1 and 20 time steps (rows) for vertical inputs 1 and 2 (cols)\n",file=outputfilepath,append=TRUE)
write.table(vertment,sep="\t",file=outputfilepath,append=TRUE,col.names=FALSE)
cat("\n",file=outputfilepath,append=TRUE)

W=c()
W[[1]]=WiEt[[20]]
W[[2]]=WEEt[[20]]
W[[3]]=WEIt[[20]]
W[[4]]=WIEt[[20]]

# Test and train vertically trained network on horizontal input
vdmulearning(W,p,r,q,s,h,50,order,horzinput,theta,hexstartendL1)

# Print layer 1 state figure at 1st, 20th, and 50th run

ent=matrix(0,3,2) # initiate entropy variable
kc=1

for (k in c(1,20,50)) {
	for (l in 1:2) {
		png(file=sprintf("horzinput_%.0f_layer1_r%.0f_fig.png",l,k))
		vdmhexplot(hexstartendL1,asEt[[k]][,l],1,1)
			
		# Compute excitatory activation for entropy
		ent[kc,l]=ent[kc,l]+sum(asEt[[k]])

		dev.off()
	} 
	kc=kc+1
}

# Calculate mean entropy
horzment=ent/169

# Print to output
cat("Mean entropy at 1, 20 and 50 time steps (rows) for horizontal inputs 1 and 2 (cols)\n",file=outputfilepath,append=TRUE)
write.table(horzment,sep="\t",file=outputfilepath,append=TRUE,col.names=FALSE)
cat("\n",file=outputfilepath,append=TRUE)


# PART 3: OBJECT CATEGORY
#-----------------------------------------------------------------

# Train with squares and crosses
W=c(19,15,5)
vdmulearning_2layer(W,p,r,q,s,h,50,order,scinput,theta,hexstartendL1,hexstartendL2,9)

# Print layer 2 state figure

ent=matrix(0,3,9) # initiate entropy variable
kc=1

for (k in c(1,20,50)) {
	for (l in 1:8) {
	
		# Print layer 1
		png(file=sprintf("scinput_%.0f_layer1_r%.0f_fig.png",l,k),width=800,height=800)
		layout(matrix(c(0,0,1,0,0,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,0,0,9,0,0),5,5,byrow=TRUE));layout.show(9)
		for (i in 1:9) {
			vdmhexplot(hexstartendL1,asEt[[k]][[1]][[i]][,l],1,1)
		} # end layer 1 retinotopic location loop
		dev.off()
		
		# Print layer 2
		png(file=sprintf("scinput_%.0f_layer2_r%.0f_fig.png",l,k))
		vdmhexplot(hexstartendL2,asEt[[k]][[2]][,l],1,1)
		dev.off()
		
	} 
	kc=kc+1
}

# Print sensitivity curves for layer 2 units at time step 50
png(file="scinput_sensitivity_layer2_fig.png",width=800,height=800)
layout(matrix(c(1:19,0),5,4,byrow=TRUE));layout.show(19)
for (u in 1:19) {
	plot(asEt[[50]][[2]][u,],type="l")
}
dev.off()
	

# END WORKSPACE COMPUTATIONS

# Close files
close(outputfile)

# Clear workspace
#rm(list=ls())
