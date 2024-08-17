#' @title Do ancestry analysis with EM algorithm
#' @description The EM algorithm could be used for estimating the Q and F matrix.
#' @usage em(g, q, f, acc, max.iter, tol, model)
#' @param g Genotype matrix with dimensions \eqn{n × p}, where n is sample size
#' and p is the number of SNPs.
#' @param q Ancestry coefficient matrix with dimensions \eqn{n × K}, where n
#' is sample size and K is the number of populations.
#' @param f Minor allele frequency matrix with dimensions \eqn{K × p},
#' where K is the number of populations and p is the number of SNPs.
#' @param acc a logical value indicating whether use accelerated EM or not.
#' @param max.iter an integer. If acc is TRUE, the number of iterations must be set.
#' @param tol Tolerance. If acc is FALSE, tol must be set. The default is 1e-4.
#' @param model Choose which model you want to use. Supervised learning or unsupervised learning.
#' @return Estimation results of q, f and the loglikelihood value of each iteration.
#' @export

# The EM algorithm
em <- function(g, q, f, acc, max.iter = 3, tol = 1e-4,
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
    p <- ncol(f)
    n <- nrow(q)
    K <- ncol(q)
    if (acc == TRUE) {
        loglikem <- numeric()
        loglikem[1] <- lfun(g, q, f)
        flongvec1 <- matrix(NA, nrow(f) * ncol(f), max.iter)
        flongvec1[, 1] <- matrix(f, nrow(f) * ncol(f), 1)
        if (model == "supervised") {
            qvec <- matrix(NA, ncol(q), max.iter)
            qvec[, 1] <- t(q[nrow(q), ])
        }
        if (model == "unsupervised") {
            qvec <- matrix(NA, nrow(q) * ncol(q), max.iter)
            qvec[, 1] <- matrix(q, nrow(q) * ncol(q), 1)
        }
        iter <- 1
        repeat{
            iter <- iter + 1
            # memory the current q and f
            q0 <- q
            f0 <- f
            if (model == "supervised") {
                q <- em.q(g, q0, q, f0, p, n)
                qvec[, iter] <- t(q[nrow(q), ])
            }
            if (model == "unsupervised") {
                for (i in 1:nrow(q0)) {
                    q[i, ] <- em.q(g, q0, q, f0, p, i)[i, ]
                }
                qvec[, iter] <- matrix(q, nrow(q) * ncol(q), 1)
            }
            # update F
            f <- .Call("updatef", n, p, K, q0, f0, g)
            f[f < lb] <- lb
            f[f > ub] <- ub
            flongvec1[, iter] <- matrix(f, nrow(f) * ncol(f), 1)
            loglikem[iter] <- lfun(g, q, f)
            if (iter == max.iter)
                break
        }
        if (model == "supervised") {
            res <- list(qvector = qvec, fvector = flongvec1,
                        loglike = loglikem, q = q, f = f)
        } else {
            res <- list(qvector = qvec, fvector = flongvec1,
                        loglike = loglikem, q = q, f = f)
        }
        return(res)
    }
    if (acc == FALSE) {
        if (is.null(tol)) {
            stop("you must set up the tolerance")
        }
        loglikem <- numeric()
        loglikem[1] <- lfun(g, q, f)
        iter <- 1
        repeat{
            iter <- iter + 1
            # memory the current q and f
            q0 <- q
            f0 <- f
            # update Q
            if (model == "supervised") {
                q <- em.q(g, q0, q, f0, p, n)
            }
            if (model == "unsupervised") {
                for (i in 1:nrow(q0)) {
                    q[i, ] <- em.q(g, q0, q, f0, p, i)[i, ]
                }
            }
            # update F
            f <- .Call("updatef", n, p, K, q0, f0, g)
            loglikem[iter] <- lfun(g, q, f)
            if (abs(loglikem[iter] - loglikem[iter - 1] < tol))
                break
        }
        return(list(f = f, loglike = loglikem[1:iter], q = q))
    }
}

#' @title EM when f is fixed
#' @description This function can be used for ancestry analysis when frequency matrix is fixed.
#' @usage fFixEm(gnew, qnew, f, acc, max.iter, tol = 1e-4, pubdata)
#' @param gnew Genotype matrix. The number of row present in gnew is 1 and the number
#' of column is the number of SNPs.
#' @param qnew Initial q used in calculation. A vector. sum(q) must be 1.
#' @param f Allele frequencies learned from the reference panels.
#' @param acc a logical value indicating whether use quasi-Newton accelerated EM or not.
#' @param max.iter an integer. If acc is TRUE, the number of iterations must be set.
#' @param tol Tolerance. If acc is FALSE, tol must be set. The default is 1e-4.
#' @param pubdata You can choose a public dataset here, E11, K13, K4, K12b, K7b, World9. You also can use other public
#' dataset which is not in this package.
#' @return Estimation results of q and the loglikelihood value of each iteration.
#' @export

fFixEm <- function(gnew, qnew, f, acc, max.iter = 3,
                   tol = 1e-4, pubdata = NULL) {
    if (is.null(acc)) {
        stop("You must set up a logical value for acc!")
    }
    p <- ncol(f)
    n <- nrow(qnew)
    if (acc == TRUE) {
        loglikem <- numeric()
        loglikem[1] <- lfun(gnew, qnew, f)
        qvec <- matrix(NA, ncol(qnew), max.iter)
        qvec[, 1] <- t(qnew[nrow(qnew), ])
        iter <- 1
        repeat {
            iter <- iter + 1
            # memory the current q and f
            q0 <- qnew
            qnew <- em.q(gnew, q0, qnew, f, p, n)
            qvec[, iter] <- t(qnew[nrow(qnew), ])
            loglikem[iter] <- lfun(gnew, qnew, f)
            if (iter == max.iter)
                break
        }
        return(list(qvector = qvec, loglike = loglikem, q = qnew))
    }
    if (acc == FALSE) {
        loglikem <- numeric()
        loglikem[1] <- lfun(gnew, qnew, f)
        iter <- 1
        repeat {
            iter <- iter + 1
            # memory the current q
            q0 <- qnew
            # update Q
            qnew <- em.q(gnew, q0, qnew, f, p, n)
            loglikem[iter] <- lfun(gnew, qnew, f)
            if (abs(loglikem[iter] - loglikem[iter - 1] < tol))
                break
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
        return(list(loglike = loglikem[1:iter], q = qnew))
    }
}
