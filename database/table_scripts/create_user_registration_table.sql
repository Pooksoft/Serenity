-- create the USER_REGISTRATION_TABLE
CREATE TABLE USER_REGISTRATION_TABLE(
    user_registration_id    INTEGER PRIMARY KEY,
    user_api_key            VARCHAR(256)    NOT NULL,
    user_email_address      TEXT            NOT NULL,
    user_organization       TEXT            NOT NULL
);

-- insert a test user into the USER_REGISTRATION_TABLE
INSERT INTO USER_REGISTRATION_TABLE (user_api_key,user_email_address,user_organization)
VALUES ('a404b4f0-e54a-44cb-9759-636338efc3c4','jvarner@pooksoft.com','Pooksoft');