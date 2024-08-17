#' @title quasi-Newton algorithm for ancestry analysis
#' @description Use quasi-Newton algorithm to accelerate EM or block relaxation.
#' @usage qn(g, q, f, tol = 1e-4, method, model)
#' @param g Genotype matrix with dimensions \eqn{n × p}, where n is sample size
#' and p is the number of SNPs.
#' @param q Ancestry coefficient matrix with dimensions \eqn{n × K}, where n
#' is sample size and K is the number of populations.
#' @param f Minor allele frequency matrix with dimensions \eqn{K × p},
#' where K is the number of populations and p is the number of SNPs.
#' @param tol Tolerance, the default value is 1e-4.
#' @param method Choose which algorithm you want to use. EM or BR.
#' @param model Choose which model you want to use. Supervised learning or unsupervised learning.
#' @return Estimation results of q, f and the loglikelihood value of each iteration.
#' @export

qn <- function(g, q, f, tol = 1e-4, method = c("EM", "BR"),
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
    em_res <- em(g, q, f, acc = T, max.iter = 5, model = model)
    if (method == "EM") {
        qvector <- em_res$qvector
        fvector <- em_res$fvector
    } else {
        br_res <- br(g, em_res$q, em_res$f, acc = T,
                     max.iter = 3, model = model)
        qvector <- br_res$qvector
        fvector <- br_res$fvector
    }
    qvector <- qvector[, (ncol(qvector) - 2):(ncol(qvector))]
    fvector <- fvector[, (ncol(fvector) - 2):(ncol(fvector))]
    iter <- 1
    loglike <- numeric()
    if (model == "supervised") {
        loglike[1] <- lfun(g, rbind(q[-nrow(q), ], qvector[, 1]),
                           matrix(fvector[, 1], nrow(f), ncol(f)))
    } else {
        loglike[1] <- lfun(g, matrix(qvector[, 1], nrow(q), ncol(q)),
                           matrix(fvector[, 1], nrow(f), ncol(f)))
    }
    q <- em_res$q
    f <- em_res$f
    repeat{
        iter <- iter + 1
        updateq <- qnq(q, qvector, model = model)
        updatef <- qnf(f, fvector)
        if (model == "supervised") {
            loglike[iter] <- lfun(g, rbind(q[-nrow(q), ], updateq), updatef)
        } else {
            loglike[iter] <- lfun(g, updateq, updatef)
        }
        if (loglike[iter] < loglike[iter - 1]) {
            if (model == "supervised") {
                q <- rbind(q[-nrow(q), ], qvector[, 3])
            } else {
                q <- matrix(qvector[, 3], nrow(q), ncol(q))
            }
            f <- matrix(fvector[, 3], nrow(f), ncol(f))
            loglike[iter] <- lfun(g, q, f)
        } else {
            if (model == "supervised") {
                q <- rbind(q[-nrow(q), ], updateq)
            } else {
                q <- updateq
            }
            f <- updatef
        }
        if (abs(lfun(g, q, f) - loglike[iter - 1]) < tol) {
            break
        } else {
            if (method == "EM") {
                fitmodel <- em(g, q, f, acc = T, max.iter = 3, model = model)
            } else {
                fitmodel <- br(g, q, f, acc = T, max.iter = 3, model = model)
            }
            qvector <- fitmodel$qvector
            fvector <- fitmodel$fvector
        }
    }
    if (model == "supervised") {
        qnew <- round(q[nrow(q), ], 5)
        res <- list(q = qnew, f = f, loglike = loglike[1:iter])
    } else {
        res <- list(q = q, f = f, loglike = loglike[1:iter])
    }
    return(res)
}


#' @title quasi-Newton when f is fixed
#' @description quasi-Newton for ancestry analysis when F is fixed
#' @usage fFixQN(gnew, qnew, f, tol, method, pubdata)
#' @param gnew Integer which length is the number of SNPs used in calculation.
#' @param qnew Initial q used in calculation. A vector. sum(q) must be 1.
#' @param f Allele frequencies learned from the reference panels.
#' @param tol Tolerance, the default value is 1e-4.
#' @param method Choose which algorithm you want to use. EM or BR.
#' @param pubdata You can choose a public dataset here, E11, K13, K4, K12b, K7b, World9. You also can use other public
#' dataset which is not in this package.
#' @return Estimation results of q and the loglikelihood value of each iteration.
#' @export

fFixQN <- function(gnew, qnew, f, tol = 1e-4,
                   method = c("EM", "BR"), pubdata = NULL) {
    if (method == "EM") {
        em_res <- fFixEm(gnew, qnew, f, acc = T, max.iter = 3)
        qvector <- em_res$qvector
    } else {
        br_res <- fFixBr(gnew, qnew, f, acc = T, max.iter = 3)
        qvector <- br_res$qvector
    }
    qvector <- qvector[, (ncol(qvector) - 2):(ncol(qvector))]
    iter <- 1
    loglike <- numeric()
    loglike[1] <- lfun(gnew, rbind(qnew[-nrow(qnew), ], qvector[, 1]), f)
    repeat{
        iter <- iter + 1
        updateq <- qnq(qnew, qvector, model = "supervised")
        loglike[iter] <- lfun(gnew, rbind(qnew[-nrow(qnew), ], updateq), f)
        if (loglike[iter] < loglike[iter - 1]) {
            qnew <- rbind(qnew[-nrow(qnew), ], qvector[, 3])
            loglike[iter] <- lfun(gnew, qnew, f)
        } else {
            qnew <- rbind(qnew[-nrow(qnew), ], updateq)
            loglike[iter] <- lfun(gnew, qnew, f)
        }
        if (abs(loglike[iter] - loglike[iter - 1]) < tol) {
            break
        } else {
            if (method == "EM") {
                fitmodel <- fFixEm(gnew, qnew, f, acc = T, max.iter = 3)
            } else {
                fitmodel <- fFixBr(gnew, qnew, f, acc = T, max.iter = 3)
            }
            qvector <- fitmodel$qvector
        }
    }
    qnew <- round(qnew, 5)
    qnew <- as.data.frame(qnew)
    rownames(qnew) <- "result"
    if(is.null(pubdata)) {
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
    return(list(q = qnew, f = f, loglike = loglike[1:iter]))
}
