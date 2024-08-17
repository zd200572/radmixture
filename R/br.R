#' @title Block Relaxation for parameters estimation
#' @description This function is also used for estimating Q and F but faster than EM.
#' @usage br(g, q, f, acc, max.iter, tol, model)
#' @param g Genotype matrix with dimensions \eqn{n × p}, where n is sample size
#' and p is the number of SNPs.
#' @param q Ancestry coefficient matrix with dimensions \eqn{n × K}, where n
#' is sample size and K is the number of populations.
#' @param f Minor allele frequency matrix with dimensions \eqn{K × p},
#' where K is the number of populations and p is the number of SNPs.
#' @param acc a logical value indicating whether use quasi-Newton accelerated BR or not.
#' @param max.iter If acc = T, max.iter must be set, the default is 3.
#' @param tol Tolerance, if acc = F, tolerance must be set, the default is 1e-4.
#' @param model Choose which model you want to use. Supervised learning or unsupervised learning.
#' @return Estimation results of q, f and the loglikelihood value of each iteration.
#' @export

# block relaxation (sequential quadratic programming)
br <- function(g, q, f, acc, max.iter = 3, tol = 1e-4,
               model = c("supervised", "unsupervised")) {
    if (!is.matrix(g)) {
        stop("g should be a matrix")
    }
    if (!is.matrix(q)) {
        stop("q should be a matrix")
    }
    if (!is.matrix(f)) {
        stop("f should be a matrix")
    }
   if (is.null(acc)) {
       stop("You must set a logical value for acc!")
   }
   if (acc == TRUE) {
        logbl <- numeric()
        logbl[1] <- lfun(g, q, f)
        flongvec1 <- matrix(NA, nrow(f) * ncol(f), max.iter)
        flongvec1[, 1] <- matrix(f, nrow(f) * ncol(f), 1)
        if (model == "supervised") {
            qlongvec1 <- matrix(NA, ncol(q), max.iter)
            qlongvec1[, 1] <- t(q[nrow(q), ])
        } else if (model == "unsupervised") {
            qlongvec1 <- matrix(NA, nrow(q) * ncol(q), max.iter)
            qlongvec1[, 1] <- matrix(q, nrow(q) * ncol(q), 1)
        }

        iter <- 1

        repeat {
            iter <- iter + 1
            f <- brf(g, q, f)
            flongvec1[, iter] <- matrix(f, nrow(f) * ncol(f), 1)
            q <- brq(g, q, f, model)
            if (model == "supervised") {
                qlongvec1[, iter] <- t(q[nrow(q), ])
            } else if (model == "unsupervised") {
                qlongvec1[, iter] <- matrix(q, nrow(q) * ncol(q), 1)
            }
            logbl[iter] <- lfun(g, q, f)
            if (iter == max.iter) {
                break
            }
        }
        res <- list(qvector = qlongvec1, fvector = flongvec1, loglike = logbl)
        return(res)
    }

    if (acc == FALSE) {
        if (is.null(tol)) {
            stop("You should set up tolerance!")
        }
        logbl <- numeric()
        logbl[1] <- lfun(g, q, f)
        iter <- 1
        repeat {
            iter <- iter + 1
            f <- brf(g, q, f)
            q <- brq(g, q, f, model)
            logbl[iter] <- lfun(g, q, f)
            if (abs(logbl[iter] - logbl[iter - 1]) < tol) {
                break
            }
        }
        if (model == "supervised") {
            res <- list(q = q[nrow(q), ], f = f, loglike = logbl[1:iter])
        } else {
            res <- list(q = q, f = f, loglike = logbl[1:iter])
        }
        return(res)
    }
}


#' @title Block relaxation when f is fixed
#' @description This function can be used for ancestry analysis when frequency matrix is fixed.
#' @usage fFixBr(gnew, qnew, f, acc, max.iter, tol, pubdata)
#' @param gnew Genotype matrix. The number of row present in gnew is 1 and the number
#' of column is the number of SNPs.
#' @param qnew Initial q used in calculation. A vector. Sum(q) must be 1.
#' @param f Allele frequencies matrix learned from the reference panels.
#' @param acc a logical value indicating whether use quasi-Newton accelerated BR or not.
#' @param max.iter If acc = T, max.iter must be set, the default is 3.
#' @param tol If acc = F, tolerance must be set, the default is 1e-4.
#' @param pubdata You can choose a public dataset here, E11, K13, K4, K12b, K7b, World9. You also can use other public
#' dataset which is not in this package.
#' @return Estimation results of q and the loglikelihood value of each iteration.
#' @export

fFixBr <- function(gnew, qnew, f, acc, max.iter = 3,
                   tol = 1e-4, pubdata = NULL) {
    if (acc == TRUE) {
        logbl <- numeric()
        logbl[1] <- lfun(gnew, qnew, f)
        qvec <- matrix(NA, ncol(qnew), max.iter)
        qvec[, 1] <- t(qnew[nrow(qnew), ])
        iter <- 1
        repeat {
            iter <- iter + 1
            qnew <- brq(gnew, qnew, f, model = "supervised")
            qvec[, iter] <- t(qnew[nrow(qnew), ])
            logbl[iter] <- lfun(gnew, qnew, f)
            if (iter == max.iter) {
                break
            }
        }
        res <- list(qvector = qvec, q = qnew[nrow(qnew), ],
                    loglike = logbl[1:iter])
        return(res)
    } else {
        if (is.null(tol)) {
            stop("You should set up the tolerance!")
        }
        logbl <- numeric()
        logbl[1] <- lfun(gnew, qnew, f)
        iter <- 1
        repeat {
            iter <- iter + 1
            qnew <- brq(gnew, qnew, f, model = "supervised")
            logbl[iter] <- lfun(gnew, qnew, f)
            if (abs(logbl[iter] - logbl[iter - 1]) < tol) {
                break
            }
        }
        qnew <- round(qnew, 5)
        qnew <- as.data.frame(qnew)
        rownames(qnew) <- "result"
        if (is.null(pubdata)) {
            colnames(qnew) <- NULL
        } else if (pubdata == "E11") {
            colnames(qnew) <- c("African", "European", "India", "Malay",
                                "SouthChineseDai", "SouthwestChineseYi",
                                "EastChinese", "Japanese", "NorthChineseOroqen",
                                "Yakut", "American")
        } else if (pubdata == "K13") {
            colnames(qnew) <- c("Siberian", "Amerindian", "West_African",
                                "Palaeo_African", "Southwest_Asian",
                                "East_Asian", "Mediterranean", "Australasian",
                                "Arctic", "West_Asian", "North_European",
                                "South_Asian", "East_African")
        } else if (pubdata == "K4") {
            colnames(qnew) <- c("European", "Asian", "African", "Amerindian")
        } else if (pubdata == "K7b") {
            colnames(qnew) <- c("South_Asian", "West_Asian", "Siberian", "African",
                                "Southern", "Atlantic_Baltic", "East_Asian") 
        } else if (pubdata == "World9") {
            colnames(qnew) <- c("Amerindian", "East_Asian", "African", "Atlantic_Baltic",
                                "Australasian", "Siberian", "Caucasus_Gedrosia", "Southern",
                                "South_Asian")
        } else if (pubdata == "K12b") {
            colnames(qnew) <- c("Gedrosia", "Siberian", "Northwest_African",
                                "Southeast_Asian", "Atlantic_Med", "North_European",
                                "South_Asian", "East_African", "Southwest_Asian",
                                "East_Asian", "Caucasus", "Sub_Saharan")
        } else if (pubdata == "Dog") {
            colnames(qnew) <- c(Bourbonnais Pointing Dog','Majorca Mastiff','Alaskan Malamute','Drever','Mudi','Cretan Tracer','Kerry Blue Terrier','Black Russian Terrier','Eurasier','Old German Shepherd','Basset Fauve de Bretagne','Griffon Fauve de Bretagne','Doberman Pinscher','South Russian Ovcharka','Dachshund','Shih Tzu','Kuvasz','Gordon Setter','Vizsla','Clumber Spaniel','Belgian Laekenois','Peruvian Hairless','Bergamasco Sheepdog','Afghan Hound','Akbash','Anatolian Shepherd Dog','Canaan Dog','Central Asian Shepherd Dog','Caucasian Ovcharka','Kangal','Saluki','Sarabi','Tazi','Turkish Mastiff','Zagar','Cane Corso','Swedish White Elkhound','Staffordshire Bull Terrier','American Staffordshire Terrier','Japanese Chin','Ariege Pointer','Rhodesian Ridgeback','Belgian Sheepdog','Groenendael','Silky Terrier','Pyrenean Shepherd','German Giant Spitz','Sealyham Terrier','Australian Cattle Dog','Australian Kelpie','Norrbottenspitz','Schipperke','Chesapeake Bay Retriever','Slovensky Kopov','Glen of Imaal Terrier','Lancashire Heeler','Boerboel','American Water Spaniel','Olde English Bulldogge','German Wirehaired Pointer','Blue Gascony Griffon','Briquet Griffon Vendeen','Grand Griffon Vendeen','Petit Basset Griffon Vendeen','Portuguese Water Dog','Auvergne Pointer','Boykin Spaniel','Weimaraner','American Foxhound','Kai Ken','Large Munsterlander','Irish Terrier','Kishu','Shiba Inu','Shikoku','Russian Tsvetnaya Bolonka','Magyar Agar','Volpino Italiano','Grand Basset Griffon Vendeen','Sarplaninac','Bluetick Coonhound','Black and Tan Coonhound','Mountain Cur','Plott Hound','Redbone Coonhound','Treeing Walker Coonhound','German Shorthaired Pointer','Bohemian Shepherd','Bouvier des Ardennes','East-European shepherd','German Shepherd Dog','Shiloh Shepherd','Tamaskan','Galgo Espanol','Havanese','Dandie Dinmont Terrier','Kars','Australian Terrier','Bolognese','Swedish Vallhund','Dogo Canario','Estrela Mountain Dog','Hanoverian Scenthound','Standard Schnauzer','German Spitz Mittel','Posavac Hound','Appenzeller Sennenhund','Cairn Terrier','Polish Tatra Sheepdog','English Shepherd','English Bulldog','Continental Bulldog','Alpine Dachsbracke','Harrier','Irish Water Spaniel','Greenland Dog','Brazilian Terrier','Swedish Lapphund','Icelandic Sheepdog','Entlebucher Mountain Dog','Hortaya Borzaya','Polish Greyhound','Flat-Coated Retriever','French Bulldog','Broholmer','Catalburun','English Setter','Leonberger','English Toy Terrier','Schapendoes','German Spitz','Ibizan Hound','Spanish Mastiff','German Pinscher','Pug','Briard','Dogo Argentino','Scottish Deerhound','Sloughi','Brittany','Finnish Hound','English Cocker Spaniel','Catalan Sheepdog','Bruno Jura Hound','Bullmastiff','Norfolk Terrier','Billy','Great Anglo-French Tricolour Hound','Norwegian Buhund','Miniature Poodle','Poodle','Toy Poodle','American Hairless Terrier','Small Munsterlander','Airedale Terrier','Ariegeois','Great Anglo-French White and Orange Hound','Grand Bleu de Gascogne','Petit Bleu de Gascogne','Porcelaine','Basset Artesien Normand','Bichon Frise','Chinese Shar-Pei','White Swiss Shepherd Dog','Turkish Zerdava','Dalmatian','Shetland Sheepdog','Bavarian Mountain Scent Hound','Spinone Italiano','Portuguese Podengo','Segugio Italiano','Pharaoh Hound','Bracco Italiano','Saint-Usuge Spaniel','Whippet','Manchester Terrier (Standard)','Field Spaniel','Pumi','Chinook','Miniature Pinscher','Prague Ratter','Sussex Spaniel','Neapolitan Mastiff','Cesky Terrier','Blue Picardy Spaniel','Picardy Spaniel','Hokkaido','Portuguese Podengo Pequeno','Cirneco dell_Etna','Labrador Retriever','Portuguese Sheepdog','Japanese Spitz','Saarloos Wolfdog','Pekingese','English Springer Spaniel','Czechoslovakian Wolfdog','Slovak Cuvac','Hallefors Elkhound','Norwegian Elkhound','Dogue de Bordeaux','Tibetan Terrier','Welsh Springer Spaniel','German Spitz Klein','Pomeranian','Tosa Inu','Pont-Audemer spaniel','Cesky Fousek','Newfoundland','Chow Chow','Norwich Terrier','Dutch Partridge Dog','Ratonero Bodeguero Andaluz','Spanish Water Dog','Xoloitzcuintli','Samoyed','English Toy Spaniel','Kromfohrlander','Gotland Hound','Hamiltonstovare','Puli','Mastiff','Hovawart','Norwegian Lundehund','Lowchen','Polish Lowland Sheepdog','Swiss Hound','Formosan Mountain Dog','Kintamani','Thai Ridgeback','German Hunting Terrier','Russian Toy','Stabyhoun','Greater Swiss Mountain Dog','Beauceron','Biewer Terrier','Yorkshire Terrier','Keeshond','Finnish Spitz','Wirehaired Pointing Griffon','Miniature Schnauzer','Giant Schnauzer','Azawakh','Chinese Crested','Miniature Bull Terrier','Unknown','Wire Fox Terrier','Welsh Terrier','Saint Bernard','Boston Terrier','Saint Miguel Cattle Dog','Picardy Shepherd','Belgian Tervuren','Belgian Malinois','American Akita','Bernese Mountain Dog','Cocker Spaniel','English Pointer','Japanese Akita','Lhasa Apso','Slovakian Wirehaired Pointer','Swedish Elkhound','Braques Francais','Portuguese Pointer','Basset Hound','Russian Hound','Danish-Swedish Farmdog','Artois Hound','Bloodhound','Hellenic Hound','Bouvier des Flandres','Affenpinscher','Irish Setter','Irish Red and White Setter','Pudelpointer','Tibetan Spaniel','Collie','Skye Terrier','American Eskimo Dog','Dutch Shepherd','Otterhound','Griffon Nivernais','French Spaniel','Pyrenean Mastiff','Chihuahua','Estonian Hound','Small Swiss Hound','Smooth Fox Terrier','English Foxhound','Elo','Rat Terrier','Teddy Roosevelt Terrier','Toy Fox Terrier','Maremma Sheepdog','Coton de Tulear','Maltese','German Hound','Brussels Griffon','Brussels Griffon (smooth)','Old English Sheepdog','Catahoula Leopard Dog','Nova Scotia Duck Tolling Retriever','Lakeland Terrier','Barbet','Parson Russell Terrier','Patterdale Terrier','Continental Toy Spaniel (Papillon)','Continental Toy Spaniel (Phalene)','American Bulldog','German Spaniel','Golden Retriever','Curly-Coated Retriever','Great Pyrenees','Landseer','Bedlington Terrier','Australian Shepherd','Boxer','Cardigan Welsh Corgi','Pembroke Welsh Corgi','Bull Terrier','Croatian Sheepdog','Silken Windhound','Cavalier King Charles Spaniel','Lagotto Romagnolo')
        }
        res <- list(q = qnew, loglike = logbl[1:iter])
        return(res)
    }
}
